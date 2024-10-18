function Data = processPAR_fit(Data,ProfileInfo,defaultPars,dataType)
% PROCESSPAR_fit is a subfunction of processPAR.m fitting a functional fit to the processed radiometric data and extending
% observations to the surface (where observations were removed during previous processing steps due to saturation). Processed
% data are appended to the existing "Data" structure.
%
% INPUT ARGUMENTS
% Data - Seal tag/BGC-Argo raw data, processed data and metadata
%   structure, created in loadData.m
% ProfileInfo - Profile specific information/metadata (general, PAR, IRR490, FLUO)
%   structure, created in loadData.m
% defaultPars - default processing parameters
%   structure, created in defaultPars.m
% dataType - Type of radiometric data
%   string, "PAR" or "IRR490"
%
% OUTPUT
% Data - Seal tag/BGC-Argo raw data, processed data and metadata
%   structure, processed data stored after each processing step

%% Shallowest and deepest available observation
firstObs = find_ndim(isfinite(Data.Processed.(dataType).log.RegDrkSat),1,'first')';
lastObs = find_ndim(isfinite(Data.Processed.(dataType).log.RegDrkSat),1,'last')';

%% Functional fit: full profile (where observations are available)
% Select profiles for fit computation
i_profiles = find(...
    ProfileInfo.General.Processed &...              % Passed processing checks in loadData
    ProfileInfo.General.Daytime &...                % Daytime profile
    sum(~isnan(Data.Processed.(dataType).log.RegDrkSat))' >= 3)';  % At least three values
% i_profiles = find(...
%     ProfileInfo.General.Processed &...              % Passed processing checks in loadData
%     ProfileInfo.General.Daytime &...                % Daytime profile
%     ProfileInfo.(dataType).PAR_QC == 1 &...
%     all(isfinite(Data.Processed.(dataType).log.RegDrkSat(21:24,:)))' & ...
%     sum(~isnan(Data.Processed.(dataType).log.RegDrkSat))' >= 3)';  % At least three values

% Fit parameters
% Initialise bspline fit coefficients
fdCoeffs = zeros(defaultPars.LFM.nBasis,1);
% Set up functional parameter object
lfdObj = 1; % original value in script = 2; penalize curvature of acceleration

% Compute fit
for iP = i_profiles
    % Dark/saturation corrected PAR data to compute fit over
    PAR_fit = Data.Processed.(dataType).log.RegDrkSat(:,iP);
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
    lastwarn('')
    [~, ~, fdFit_PAR, ~, ~, ~, ~] = smooth_monotone(fitInterval,PAR_fit,fdPar_iP);
    % Catch warning that PAR data is badly scaled for monotone fit, do not proceed to calculate fitted data in that case
    [~, warnID] = lastwarn;
    if ~strcmp(warnID, 'MATLAB:nearlySingularMatrix')
        % Evaluate fit and store in platform_processed PAR structure
        Data.Processed.(dataType).log.RegDrkSatFitAll(fitInterval,iP) = eval_fd(fdFit_PAR,fitInterval);

        % Calculate Kd as the derivative function of the PAR fit
        fdKd_iP = deriv_fd(fdFit_PAR);
        Data.Processed.(dataType).Kd.FitAll(fitInterval,iP) = -eval_fd(fdKd_iP,fitInterval);
    else
        % If warning was issued, ignore profile in the following processing steps
        i_profiles(i_profiles==iP) = [];
    end
end

% Compute linear PAR array
Data.Processed.(dataType).lin.RegDrkSatFitAll = exp(Data.Processed.(dataType).log.RegDrkSatFitAll);

%% Extend fitted data to the surface
% Linearly extrapolate log-transformed light, assuming continuous exponential increase of light to the surface
Data.Processed.(dataType).log.RegDrkSatFitSurf = Data.Processed.(dataType).log.RegDrkSatFitAll;

for iP = i_profiles
    % Depth range over which to calculate the log-transformed increase of light: 
    % >> first available observation to MLD (at least 5 observations in the ML required, else no extrapolation performed)
    lightExtrapDepth = abs(round(ProfileInfo.General.MLD(iP)));
    if isfinite(lightExtrapDepth) && lightExtrapDepth >= firstObs(iP) + 4
        % Linear model through light vs depth over the defined depth range
        PAR_lmfit = fitlm( ...
            -1 : -1 : -lightExtrapDepth, ...
            Data.Processed.(dataType).log.RegDrkSatFitAll(1:lightExtrapDepth,iP));

        Data.Processed.(dataType).log.RegDrkSatFitSurf(1,iP) = Data.Processed.(dataType).log.RegDrkSatFitSurf(firstObs(iP),iP) + PAR_lmfit.Coefficients.Estimate(2)*(firstObs(iP)-1);
        Data.Processed.(dataType).log.RegDrkSatFitSurf(:,iP) = fillmissing(Data.Processed.(dataType).log.RegDrkSatFitSurf(:,iP),'linear','EndValues','none');
    end
end

% Compute linear PAR array
Data.Processed.(dataType).lin.RegDrkSatFitSurf = exp(Data.Processed.(dataType).log.RegDrkSatFitSurf);

%% Fit another bspline to the surface-extended data for smoothing
Data.Processed.(dataType).log.RegDrkSatFitSurfSmooth = Data.Processed.(dataType).log.RegDrkSatFitSurf;

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
    PAR_fit = Data.Processed.(dataType).log.RegDrkSatFitSurf(:,iP);
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
    Data.Processed.(dataType).log.RegDrkSatFitSurfSmooth(fitInterval,iP) = eval_fd(fdFit_PAR,fitInterval);

    % % Calculate Kd as the derivative function of the PAR fit
    % fdKd_iP = deriv_fd(fdFit_PAR);
    % Data.Processed.(dataType).Kd.FitAll(fitInterval,iP) = -eval_fd(fdKd_iP,fitInterval);
end

% Compute linear PAR array
Data.Processed.(dataType).lin.RegDrkSatFitSurfSmooth = exp(Data.Processed.(dataType).log.RegDrkSatFitSurfSmooth);

%% Functional fit (predefined depth interval)
% Select profiles for fit computation
lumToFit_bnd = Data.Processed.(dataType).log.RegDrkSat;
i_profiles = find(...
    (firstObs <= abs(defaultPars.upperDepthBound)...    % Finite data within bounds
    & lastObs >= abs(defaultPars.lowerDepthBound)) &...
    ProfileInfo.General.Processed &...                          % Passed processing checks in loadData
    ProfileInfo.General.Daytime);                               % Daytime profile

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
    Data.Processed.(dataType).log.RegDrkSatFitBnd(fitInterval,i_profiles) = eval_fd(fdFit_PAR,fitInterval);

    % Calculate Kd as the derivative function of the PAR fit
    fdKd = deriv_fd(fdFit_PAR);
    Data.Processed.(dataType).Kd.FitBnd(fitInterval,i_profiles) = -eval_fd(fdKd,fitInterval);

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
    Data.Processed.CHL_LFM(fitInterval,i_profiles) = eval_fd(newObj_fd,fitInterval);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREDICTION POST-PROCESSING
    % remove negative values
    Data.Processed.CHL_LFM(Data.Processed.CHL_LFM < 0) = NaN;
    % fill upper part of profile
    % CHLAupper = repmat(platform_processed.CHL_LFM(topBound_lfm,:),topBound_lfm - 1,1);
    % CHLA_LFM(1:topBound_lfm - 1,:) = CHLAupper;
end

% Compute linear PAR array
Data.Processed.(dataType).lin.RegDrkSatFitBnd = exp(Data.Processed.(dataType).log.RegDrkSatFitBnd);