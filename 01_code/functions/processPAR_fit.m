function tagProcessed = processPAR_fit(tagProcessed,ProfileInfo,ProfileInfo_PAR,defaultPars)

%% Shallowest and deepest available observation
firstObs = find_ndim(isfinite(tagProcessed.PAR.log.RegDrkSat),1,'first')';
lastObs = find_ndim(isfinite(tagProcessed.PAR.log.RegDrkSat),1,'last')';

%% Functional fit: full profile (where observations are available)
% Select profiles for fit computation
i_profiles = find(...
    ProfileInfo.Processed &...              % Passed processing checks in loadData
    ProfileInfo.Daytime &...                % Daytime profile
    sum(~isnan(tagProcessed.PAR.log.RegDrkSat))' >= 3)';  % At least three values
% i_profiles = find(...
%     ProfileInfo.Processed &...              % Passed processing checks in loadData
%     ProfileInfo.Daytime &...                % Daytime profile
%     ProfileInfo_PAR.PAR_QC == 1 &...
%     all(isfinite(tagProcessed.PAR.log.RegDrkSat(21:24,:)))' & ...
%     sum(~isnan(tagProcessed.PAR.log.RegDrkSat))' >= 3)';  % At least three values

% Fit parameters
% Initialise bspline fit coefficients
fdCoeffs = zeros(defaultPars.LFM.nBasis,1);
% Set up functional parameter object
lfdObj = 1; % original value in script = 2; penalize curvature of acceleration

% Compute fit
for iP = i_profiles
    % Dark/saturation corrected PAR data to compute fit over
    PAR_fit = tagProcessed.PAR.log.RegDrkSat(:,iP);
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % PAR_fit(1:20) = NaN;
    % firstObs(iP) = find_ndim(isfinite(PAR_fit),1,'first')';
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fitInterval = find(isfinite(PAR_fit));
    PAR_fit = PAR_fit(fitInterval);

    % Bspline functional fit
    basisObj = create_bspline_basis([fitInterval(1) fitInterval(end)],defaultPars.LFM.nBasis,defaultPars.LFM.nOrder);
    fdObj = fd(fdCoeffs, basisObj);
    lambda = defaultPars.LFM.lambda; %  smoothing parameter %10^(-0.1), Q: Why is lambda=0.03 in LFM file that's imported and used in the previous fit???
    fdPar_iP = fdPar(fdObj,lfdObj,lambda);

    % Monotone smooth fit
    [~, ~, fdFit_PAR, ~, ~, ~, ~] = smooth_monotone(fitInterval,PAR_fit,fdPar_iP);

    % Evaluate fit and store in platform_processed PAR structure
    tagProcessed.PAR.log.RegDrkSatFitAll(fitInterval,iP) = eval_fd(fdFit_PAR,fitInterval);

    % Calculate Kd as the derivative function of the PAR fit
    fdKd_iP = deriv_fd(fdFit_PAR);
    tagProcessed.PAR.Kd.FitAll(fitInterval,iP) = -eval_fd(fdKd_iP,fitInterval);
end

% Compute linear PAR array
tagProcessed.PAR.lin.RegDrkSatFitAll = exp(tagProcessed.PAR.log.RegDrkSatFitAll);

%% Extend fitted data to the surface
% Linearly extrapolate log-transformed light, assuming continuous exponential increase of light to the surface
tagProcessed.PAR.log.RegDrkSatFitSurf = tagProcessed.PAR.log.RegDrkSatFitAll;

for iP = i_profiles
    % Depth range over which to calculate the log-transformed increase of light: 
    % >> first available observation to MLD (at least 5 observations in the ML required, else no extrapolation performed)
    lightExtrapDepth = abs(round(ProfileInfo.MLD(iP)));
    if isfinite(lightExtrapDepth) && lightExtrapDepth >= firstObs(iP) + 4
        % Linear model through light vs depth over the defined depth range
        PAR_lmfit = fitlm( ...
            -1 : -1 : -lightExtrapDepth, ...
            tagProcessed.PAR.log.RegDrkSatFitAll(1:lightExtrapDepth,iP));

        tagProcessed.PAR.log.RegDrkSatFitSurf(1,iP) = tagProcessed.PAR.log.RegDrkSatFitSurf(firstObs(iP),iP) + PAR_lmfit.Coefficients.Estimate(2)*(firstObs(iP)-1);
        tagProcessed.PAR.log.RegDrkSatFitSurf(:,iP) = fillmissing(tagProcessed.PAR.log.RegDrkSatFitSurf(:,iP),'linear','EndValues','none');
    end
end

% Compute linear PAR array
tagProcessed.PAR.lin.RegDrkSatFitSurf = exp(tagProcessed.PAR.log.RegDrkSatFitSurf);

%% Fit another bspline to the surface-extended data for smoothing
tagProcessed.PAR.log.RegDrkSatFitSurfSmooth = tagProcessed.PAR.log.RegDrkSatFitSurf;

% Fit parameters
% Initialise bspline fit coefficients
fdCoeffs = zeros(defaultPars.LFM.nBasis,1);
% Set up functional parameter object
lfdObj = 1; % original value in script = 2; penalize curvature of acceleration

% Compute fit
for iP = i_profiles
    % Use fitted PAR data from the previous step and linear extrapolate values to the surface (keep NaNs at depth)
    % The linear extrapolation approximates the continuing increase of light to the surface and will be improved by a
    % second spline fit.
    PAR_fit = tagProcessed.PAR.log.RegDrkSatFitSurf(:,iP);
    fitInterval = find(isfinite(PAR_fit));
    PAR_fit = PAR_fit(fitInterval);

    % Bspline functional fit
    basisObj = create_bspline_basis([fitInterval(1) fitInterval(end)],defaultPars.LFM.nBasis,defaultPars.LFM.nOrder);
    fdObj = fd(fdCoeffs, basisObj);
    lambda = defaultPars.LFM.lambda; %  smoothing parameter %10^(-0.1), Q: Why is lambda=0.03 in LFM file that's imported and used in the previous fit???
    fdPar_iP = fdPar(fdObj,lfdObj,lambda);

    % Monotone smooth fit
    [~, ~, fdFit_PAR, ~, ~, ~, ~] = smooth_monotone(fitInterval,PAR_fit,fdPar_iP);

    % Evaluate fit and store in platform_processed PAR structure
    tagProcessed.PAR.log.RegDrkSatFitSurfSmooth(fitInterval,iP) = eval_fd(fdFit_PAR,fitInterval);

    % % Calculate Kd as the derivative function of the PAR fit
    % fdKd_iP = deriv_fd(fdFit_PAR);
    % tagProcessed.PAR.Kd.FitAll(fitInterval,iP) = -eval_fd(fdKd_iP,fitInterval);
end

% Compute linear PAR array
tagProcessed.PAR.lin.RegDrkSatFitSurfSmooth = exp(tagProcessed.PAR.log.RegDrkSatFitSurfSmooth);

%% Functional fit (predefined depth interval)
% Select profiles for fit computation
lumToFit_bnd = tagProcessed.PAR.log.RegDrkSat;
i_profiles = find(...
    (firstObs <= abs(defaultPars.upperDepthBound)...    % Finite data within bounds
    & lastObs >= abs(defaultPars.lowerDepthBound)) &...
    ProfileInfo.Processed &...                          % Passed processing checks in loadData
    ProfileInfo.Daytime);                               % Daytime profile

% Fit paramters
% Initialise bspline fit coefficients
fdCoeffs = zeros(defaultPars.LFM.nBasis,length(i_profiles));
% Set up functional parameter object
lfdObj = 1;

% proceed only if any profiles were found to have finite PAR data within the specified bounds
if ~isempty(i_profiles)
    % PAR data to compute fit over
    fitInterval = (abs(defaultPars.upperDepthBound) : abs(defaultPars.lowerDepthBound))';
    PAR_fit = lumToFit_bnd(fitInterval,i_profiles);

    % Bspline functional fit
    basisObj = create_bspline_basis([fitInterval(1) fitInterval(end)],defaultPars.LFM.nBasis,defaultPars.LFM.nOrder);
    fdObj = fd(fdCoeffs, basisObj);
    lambda = 0.03; %  smoothing parameter %10^(-0.1), Q: Why is lambda=0.03 in LFM file that's imported and used in the previous fit???
    fdPar_iP = fdPar(fdObj,lfdObj,lambda);

    % Monotone smooth fit
    [~, ~, fdFit_PAR, ~, ~, ~, ~] = smooth_monotone(fitInterval,PAR_fit,fdPar_iP);

    % Evaluate fit and store in PAR
    tagProcessed.PAR.log.RegDrkSatFitBnd(fitInterval,i_profiles) = eval_fd(fdFit_PAR,fitInterval);

    % Calculate Kd as the derivative function of the PAR fit
    fdKd = deriv_fd(fdFit_PAR);
    tagProcessed.PAR.Kd.FitBnd(fitInterval,i_profiles) = -eval_fd(fdKd,fitInterval);

    %% LFM CHL PREDICTION FROM KD
    % recover Blind coeffs
    KdFDcoefs(:,i_profiles) = getcoef(fdKd);
    lumBlindForPred_fd = fd(KdFDcoefs(:,i_profiles),basisObj);

    % compute Blind prediction
    % CENTER THE PREDICTOR DATA
    xcBlind_fd = center(lumBlindForPred_fd);
    % COMPUTE PREDICTION
    yhatpenBlind_fdCoefs_temp = defaultPars.LFM.LFMcore.Bpen*getcoef(xcBlind_fd);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE CHLA PREDICTION
    % eval CHLAFITBlind / CHLAHATpBlind
    newObj_dfCoefs = yhatpenBlind_fdCoefs_temp + defaultPars.LFM.LFMcore.meanFittedObsFDcoefs;
    newObj_fd = fd(newObj_dfCoefs, basisObj);
    tagProcessed.CHL_LFM(fitInterval,i_profiles) = eval_fd(newObj_fd,fitInterval);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREDICTION POST-PROCESSING
    % remove negative values
    tagProcessed.CHL_LFM(tagProcessed.CHL_LFM < 0) = NaN;
    % fill upper part of profile
    % CHLAupper = repmat(platform_processed.CHL_LFM(topBound_lfm,:),topBound_lfm - 1,1);
    % CHLA_LFM(1:topBound_lfm - 1,:) = CHLAupper;
end

% Compute linear PAR array
tagProcessed.PAR.lin.RegDrkSatFitBnd = exp(tagProcessed.PAR.log.RegDrkSatFitBnd);