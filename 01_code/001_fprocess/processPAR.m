function [platform_processed,parData] = processPAR(platform_metadata,platform_processed,genData,defVals)
% PROCESSPAR
%
% INPUT ARGUMENTS
%
% OUTPUT

%% CMD message: start
fprintf('Processing <strong>PAR</strong> data...');

%% Create parData info table
var_names = {...
    'Profile', ...              % profile number
    'noData',...                % no data in profile
    'surfaceValue', ...         % subsurface light value (first non NaN value of light)
    'surfaceValueFit', ...      % subsurface light value of fitted data
    'darkDepth', ...            % index/depth of first dark PAR value (index, depth)
    'darkValue', ...            % dark PAR value (median of PAR below dark depth)
    'SaturationDepth', ...      % depth above which light sensor maxes out
    'SaturationValue', ...      % PAR value at which light sensor maxes out
    'quenchDepth', ...          % depth of quenching threshold 15 umol m-2 s-1 (Xing et al. 2018)
    'Zeu', ...                  % euphotic depth (Morel and Berthon, 1989)
    'ZeuInterp', ...            % interpolated/smoothed euphotic depths
    'FirstOD', ...              % First optical depth
    'meanKd', ...               % mean attenuation coefficient Kd (calculated where LIGHT is non DARK non SAT)
    'attSlope', ...             % attenuation slope
    'attSlopeBound', ...        % attenuation slope between pre-defined upper and lower depth bound
    'predIntChla' ...           % prediction of integrated CHLA from light attenuation slope linear regression
    }';
parData = array2table(NaN(platform_metadata.nProfs, numel(var_names)),'VariableNames',var_names);

%% Populate initial information
% Profile #
parData.Profile = (1 : platform_metadata.nProfs)';

% Find profiles without usable data (all NaNs)
parData.noData = all(isnan(platform_processed.PAR.log_Reg))';

% Shallowest and deepest available observation
firstObs = find_ndim(isfinite(platform_processed.PAR.log_Reg),1,'first')';
lastObs = find_ndim(isfinite(platform_processed.PAR.log_Reg),1,'last')';

% subsurface light value (first non-NaN value)
nonZeroIndices = find(firstObs(:,1) ~= 0);
parData.surfaceValue(nonZeroIndices,1) = arrayfun(@(x) platform_processed.PAR.lin_Reg(firstObs(x,1), x), nonZeroIndices);

%% Dark depth/value & attenuation slope calculation
PAR_lin_RegDrk = platform_processed.PAR.lin_Reg;
PAR_log_RegDrk = platform_processed.PAR.log_Reg;

% Define daily profile subsets for dark correction calculation
% Subsets comprise the n deepest profiles measured during each deployment day (n is defined in defaults)
[uniqueVal,uniqueInd] = unique(genData.DeployDay);

% Get indices of the deepest profiles of each deployment day (indices start at 1 for each day)
[~,i_deepestPerDay] = arrayfun(...
    @(a) mink(lastObs(genData.DeployDay == a),defVals.PAR.nProfilesPerDay_dark),...
    uniqueVal,'UniformOutput',false);

% Add deployment day to indices to get accurate list of indices for the whole deployment
i_deepestPerDay = arrayfun(...
    @(a) i_deepestPerDay{a} + uniqueInd(a) - 1,...
    1:numel(uniqueVal),'UniformOutput',false);

% Sort and convert to logical indices of profiles to be used in the dark correction calculations
i_deepestPerDay = sort(vertcat(i_deepestPerDay{1:end}));
i_deepestPerDay = ismember(1:platform_metadata.nProfs,i_deepestPerDay);

% Loop through profiles selected for processing
for iP = find(i_deepestPerDay)
    %% 1) Detect and remove constant part of the PAR profile to avoid normality test on constant data series
    % See this example:
    % [h_temp,p_temp] = lillietest(ones(1,defVals.depthInterpGrid(end)),'Alpha',0.01,'Distribution','normal');
    % returns h_temp = 1, i.e. rejection of the null hypothesis (= distribution is normal)
    
    % iteratively move down the water column and check the range of PAR values between depth iZ and the deepest PAR value
    for iZ = firstObs(iP) : lastObs(iP) - 1

        iZ_PAR_range = range(PAR_lin_RegDrk(iZ:lastObs(iP),iP));

        if iZ_PAR_range <= defVals.PAR.delta_constantPAR
            Z_const = iZ + 1;           % depth of first value belonging to constant signal
            PAR_lin_RegDrk(Z_const:end,iP) = NaN;
            PAR_log_RegDrk(Z_const:end,iP) = NaN;
            lastObs(iP) = Z_const-1;    % Update last available observation depth

            break   % end for-loop when PAR range falls below threshold value where PAR is considered constant
        end

    end

    %% 2) Detect PAR dark value (Organelli et al., 2016)
    % initialize values
    Z_dark = NaN;       % dark PAR depth
    PAR_dark = NaN;     % dark PAR value
    
    % Sampled vector X in Lilliefors test must have at least 4 valid observations
    if lastObs(iP) - firstObs(iP) + 1 >= 4

        % iteratively move down the water column and check normality of distribution of PAR values from depth i and the deepest PAR value
        for iZ = firstObs(iP) : lastObs(iP) - 3
            warning('off','stats:lillietest:OutOfRangePLow')    % Warning: P is less than the smallest tabulated value, returning 0.001.
            warning('off','stats:lillietest:OutOfRangePHigh')
            
            % Lilliefor normality test
            [h_lillie,~] = lillietest(PAR_lin_RegDrk(iZ:lastObs(iP),iP),'Alpha',0.01,'Distribution','normal');
            
            if h_lillie == 0
                % If Lilliefor test returns h = 0 (i.e. distribution of PAR values from iZ to deepest value is normal)
                % then iZ = dark depth and get dark PAR value at dark depth
                Z_dark = iZ;
                PAR_dark = median(PAR_lin_RegDrk(Z_dark:lastObs(iP),iP));

                break % end for-loop
            end
        end

    end
   
    % Write to parData table
    parData.darkDepth(iP) = -Z_dark;     % first dark depth (index, depth)
    parData.darkValue(iP) = PAR_dark;	% dark PAR value

    %% 3) Attenuation slope (based on logarithmic PAR values)
    % Slope calculated over non-dark section of PAR profile
    if isfinite(Z_dark) && firstObs(iP) > 0
        atten_slope = ...
            (PAR_log_RegDrk(Z_dark-1,iP) - PAR_log_RegDrk(firstObs(iP),iP)) /...
            (defVals.depthInterpGrid(Z_dark-1) - defVals.depthInterpGrid(firstObs(iP)));

        % Slope calculated between defined upper and lower depth bound
        % (only if dark depth is deeper than lower bound)
        if defVals.lowerDepthBound > parData.darkDepth(iP)
            atten_slope_PART = ...
                (PAR_log_RegDrk(abs(defVals.lowerDepthBound),iP) - PAR_log_RegDrk(firstObs(iP),iP)) /...
                (platform_processed.DEPTH(abs(defVals.lowerDepthBound),iP) - platform_processed.DEPTH(firstObs(iP),iP));
        else
            atten_slope_PART = NaN;
        end

        % write parData table
        parData.attSlope(iP)        = atten_slope;	% slope of attenuation (calculated where PAR is non DARK)
        parData.attSlopeBound(iP)    = atten_slope_PART; % slope of attenuation calculated on [upperDepthBound:lowerDepthBound] interval
    end

end   

% extend computed values of selected profiles to remaining profiles
parData.darkDepth = fillmissing(parData.darkDepth,'nearest');
parData.darkValue = fillmissing(parData.darkValue,'nearest');

%% Dark correction
% Set PAR below dark depth NaN
darkMask = false(size(PAR_lin_RegDrk));
for iP = 1 : platform_metadata.nProfs
    darkMask(abs(parData.darkDepth(iP,1)):end,iP) = true;
end
PAR_lin_RegDrk(darkMask) = NaN;

% Offset all non-dark PAR values by dark value (profile-specific)
darkValues = repmat(parData.darkValue.',size(PAR_lin_RegDrk,1),1);
PAR_lin_RegDrk = PAR_lin_RegDrk - darkValues;

% Remove negative linear PAR values and interpolate linearly
PAR_lin_RegDrk(PAR_lin_RegDrk<=0) = NaN;
PAR_lin_RegDrk = fillmissing(PAR_lin_RegDrk,'linear',1,'EndValues','none');

% Write to platform_processed structure
platform_processed.PAR.lin_RegDrk = PAR_lin_RegDrk;
platform_processed.PAR.log_RegDrk = log(PAR_lin_RegDrk);

% Update last non-NaN parData table column
lastObs = find_ndim(isfinite(platform_processed.PAR.log_RegDrk),1,'last')';

%% NOT SURE IF ATTENUATION SLOPE SHOULD BE CALCULATED BEFORE OR AFTER DARK CORRECTION
% ???????????????????????????????????????
% ???????????????????????????????????????
% ???????????????????????????????????????
% ???????????????????????????????????????
% ???????????????????????????????????????
% ???????????????????????????????????????

%% Saturation depth
PAR_lin_RegDrkSat = platform_processed.PAR.lin_RegDrk;

for iP = 1 : platform_metadata.nProfs
    % derivative of log dark-corrected PAR
    diffPARlog = diff(PAR_log_RegDrk(:,iP));
    % mean of derivative
    meandiffPARlog = mean(diffPARlog,'omitnan');
    
    % check saturation depth only if profile starts with descrease>mean derivative
    % else no saturation depth and saturation value is defined to maximum PAR value measured
    if isfinite(meandiffPARlog) && diffPARlog(firstObs(iP)) > meandiffPARlog
        Z_saturation = firstObs(iP) + 1;
        while diffPARlog(Z_saturation) >= meandiffPARlog
            Z_saturation = Z_saturation + 1;
        end
        % Remove saturated values
        PAR_lin_RegDrkSat(1 : Z_saturation-1,iP) = NaN;
        
        % Write saturation depth and PAR value to parData table
        parData.SaturationDepth(iP) = -Z_saturation;
        parData.SaturationValue(iP) = PAR_lin_RegDrkSat(Z_saturation,iP);
    end
end

% Write to platform_processed structure
platform_processed.PAR.lin_RegDrkSat = PAR_lin_RegDrkSat;
platform_processed.PAR.log_RegDrkSat = log(PAR_lin_RegDrkSat);

% Update first non-NaN parData table column
firstObs = find_ndim(isfinite(platform_processed.PAR.lin_RegDrkSat),1,'first')';

%% Functional fit (full profile)
platform_processed.PAR.log_RegDrkSatFitAll = NaN(size(platform_processed.PAR.log_RegDrkSat));
platform_processed.PAR.Kd_all = NaN(size(platform_processed.PAR.log_RegDrkSat));
platform_processed.CHL_LFM.Full = NaN(size(platform_processed.PAR.log_RegDrkSat));

% Select profiles for fit computation
lumToFit_all = platform_processed.PAR.log_RegDrkSat;
i_profiles = find(...
    genData.Processed &...              % Passed processing checks in loadData
    genData.Daytime &...                % Daytime profile
    sum(~isnan(platform_processed.PAR.log_RegDrkSat))' >= 3)';  % At least three values

% Fit parameters
% Initialise bspline fit coefficients
fdCoeffs = zeros(defVals.LFM.nBasis,1);
% Set up functional parameter object
lfdObj = 1; % original value in script = 2; penalize curvature of acceleration

% Compute fit
for iP = i_profiles
    % PAR data to compute fit over
    PAR_fit = lumToFit_all(:,iP);
    fitInterval = find(isfinite(PAR_fit));
    PAR_fit = PAR_fit(fitInterval);

    % Bspline functional fit
    basisObj = create_bspline_basis([fitInterval(1) fitInterval(end)],defVals.LFM.nBasis,defVals.LFM.nOrder);
    fdObj = fd(fdCoeffs, basisObj);
    lambda = defVals.LFM.lambda; %  smoothing parameter %10^(-0.1), Q: Why is lambda=0.03 in LFM file that's imported and used in the previous fit???
    fdPar_iP = fdPar(fdObj,lfdObj,lambda);

    % Monotone smooth fit and get fit coefficients
    [~, ~, fdFit_PAR, ~, ~, ~, ~] = smooth_monotone(fitInterval,PAR_fit,fdPar_iP);
    
    % Evaluate fit and store in platform_processed PAR structure
    platform_processed.PAR.log_RegDrkSatFitAll(fitInterval,iP) = eval_fd(fdFit_PAR,fitInterval);

    % Calculate Kd as the derivative function of the PAR fit
    fdKd_iP = deriv_fd(fdFit_PAR);
    platform_processed.PAR.Kd_all(fitInterval,iP) = -eval_fd(fdKd_iP,fitInterval);

    %% LFM CHL PREDICTION FROM KD
    % recover Blind coeffs
    KdFDcoefs = getcoef(fdKd_iP);
    lumBlindForPred_fd = fd(KdFDcoefs,basisObj);

    % compute Blind prediction
    % CENTER THE PREDICTOR DATA
    xcBlind_fd = center(lumBlindForPred_fd);
    % COMPUTE PREDICTION
    yhatpenBlind_fdCoefs_temp = defVals.LFM.LFMcore.Bpen*getcoef(xcBlind_fd);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE CHLA PREDICTION
    % eval CHLAFITBlind / CHLAHATpBlind
    newObj_dfCoefs = yhatpenBlind_fdCoefs_temp + defVals.LFM.LFMcore.meanFittedObsFDcoefs;
    newObj_fd = fd(newObj_dfCoefs, basisObj);
    platform_processed.CHL_LFM.Full(fitInterval,iP) = eval_fd(newObj_fd,fitInterval);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREDICTION POST-PROCESSING
    % remove negative values
    platform_processed.CHL_LFM.Full(platform_processed.CHL_LFM.Full < 0) = NaN;

end

% Compute linear PAR array
platform_processed.PAR.lin_RegDrkSatFitAll = exp(platform_processed.PAR.log_RegDrkSatFitAll);

%% Functional fit (predefined depth interval)
platform_processed.PAR.log_RegDrkSatFitBnd = NaN(size(platform_processed.PAR.log_RegDrkSat));
platform_processed.PAR.Kd_bnd = NaN(size(platform_processed.PAR.log_RegDrkSat));
platform_processed.CHL_LFM.Bound = NaN(size(platform_processed.PAR.log_RegDrkSat));

% Select profiles for fit computation
lumToFit_bnd = platform_processed.PAR.log_RegDrkSat;
i_profiles = find(...
    (firstObs <= abs(defVals.upperDepthBound)...    % Finite data within bounds
    & lastObs >= abs(defVals.lowerDepthBound)) &...
    genData.Processed &...                      % Passed processing checks in loadData
    genData.Daytime);                           % Daytime profile
    
% Fit paramters
% Initialise bspline fit coefficients
fdCoeffs = zeros(defVals.LFM.nBasis,length(i_profiles));
% Set up functional parameter object
lfdObj = 1;

% proceed only if any profiles were found to have finite PAR data within the specified bounds
if ~isempty(i_profiles)
    % PAR data to compute fit over
    fitInterval = (abs(defVals.upperDepthBound) : abs(defVals.lowerDepthBound))';
    PAR_fit = lumToFit_bnd(fitInterval,i_profiles);

    % Bspline functional fit
    basisObj = create_bspline_basis([fitInterval(1) fitInterval(end)],defVals.LFM.nBasis,defVals.LFM.nOrder);
    fdObj = fd(fdCoeffs, basisObj);
    lambda = 0.08; %  smoothing parameter %10^(-0.1), Q: Why is lambda=0.03 in LFM file that's imported and used in the previous fit???
    fdPar_iP = fdPar(fdObj,lfdObj,lambda);

    % Monotone smooth fit and get fit coefficients
    [~, ~, fdFit_PAR, ~, ~, ~, ~] = smooth_monotone(fitInterval,PAR_fit,fdPar_iP);

    % Evaluate fit and store in PAR
    platform_processed.PAR.log_RegDrkSatFitBnd(fitInterval,i_profiles) = eval_fd(fdFit_PAR,fitInterval);

    % Calculate Kd as the derivative function of the PAR fit
    fdKd = deriv_fd(fdFit_PAR);
    platform_processed.PAR.Kd_bnd(fitInterval,i_profiles) = -eval_fd(fdKd,fitInterval);

    %% LFM CHL PREDICTION FROM KD
    % recover Blind coeffs
    KdFDcoefs(:,i_profiles) = getcoef(fdKd);
    lumBlindForPred_fd = fd(KdFDcoefs(:,i_profiles),basisObj);

    % compute Blind prediction
    % CENTER THE PREDICTOR DATA
    xcBlind_fd = center(lumBlindForPred_fd);
    % COMPUTE PREDICTION
    yhatpenBlind_fdCoefs_temp = defVals.LFM.LFMcore.Bpen*getcoef(xcBlind_fd);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE CHLA PREDICTION
    % eval CHLAFITBlind / CHLAHATpBlind
    newObj_dfCoefs = yhatpenBlind_fdCoefs_temp + defVals.LFM.LFMcore.meanFittedObsFDcoefs;
    newObj_fd = fd(newObj_dfCoefs, basisObj);
    platform_processed.CHL_LFM.Bound(fitInterval,i_profiles) = eval_fd(newObj_fd,fitInterval);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREDICTION POST-PROCESSING
    % remove negative values
    platform_processed.CHL_LFM.Bound(platform_processed.CHL_LFM.Bound < 0) = NaN;
    % fill upper part of profile
    % CHLAupper = repmat(platform_processed.CHL_LFM.Bound(topBound_lfm,:),topBound_lfm - 1,1);
    % CHLA_LFM(1:topBound_lfm - 1,:) = CHLAupper;
end

% Compute linear PAR array
platform_processed.PAR.lin_RegDrkSatFitBnd = exp(platform_processed.PAR.log_RegDrkSatFitBnd);

%% Final things
% - New subsurface PAR value based on fitted data
% - Euphotic depth (1% of the surface PAR depth, after Morel and Berthon, 1989)
iSurfVal_fit = find_ndim(isfinite(platform_processed.PAR.lin_RegDrkSatFitAll),1,'first');
finiteIndices = find(iSurfVal_fit>0);
for iP = finiteIndices
    parData.surfaceValueFit(iP,1) = platform_processed.PAR.lin_RegDrkSatFitAll(iSurfVal_fit(iP),iP);
    parData.Zeu(iP,1) = -find(platform_processed.PAR.lin_RegDrkSatFitAll(:,iP) / parData.surfaceValueFit(iP) < 0.01,1,'first');    
end
% fill missing and smooth Zeu data (window size: ~n profiles per day)
parData.ZeuInterp = fillmissing(parData.Zeu,'linear','EndValues','none');
parData.ZeuInterp = movmedian(parData.ZeuInterp,ceil(platform_metadata.nProfs / max(genData.DeployDay)),'omitnan','Endpoints','shrink');

% penetration depth (first optical depth)
parData.FirstOD = parData.ZeuInterp / 4.6;

% depth of quenching threshold 15 umol m-2 s-1 (Xing et al. 2018)
iZ_quench = find_ndim(platform_processed.PAR.lin_RegDrkSatFitAll > 15, 1, 'last')';
iZ_quench(iZ_quench == 0) = NaN;
parData.quenchDepth = -iZ_quench;

% (ML-)mean Kd
parData.meanKd = mean(platform_processed.PAR.Kd_all,1,'omitnan')';
for iP = 1 : platform_metadata.nProfs
    iZ_MLD = abs(round(genData.MLD(iP)));
    parData.meanKdML(iP) = mean(platform_processed.PAR.Kd_all(1:iZ_MLD,iP),1,'omitnan');
end

%% CMD message: done
fprintf('\b\b \x2713\n')

end