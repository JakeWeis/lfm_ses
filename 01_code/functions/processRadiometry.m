function [Data,ProfileInfo] = processRadiometry(Data,ProfileInfo,defaultPars,dataType)
% PROCESSRADIOMETRY processes radiometric seal tag/BGC-Argo data. Processed data and radiometric profile information are appended to
% the existing structures "Data" and "ProfileInfo".
%
% Radiometric profile information:
% 'Profile'              profile number
% 'noData'               no data in profile
% 'surfaceValue'         subsurface light value (first non NaN value of light)
% 'surfaceValueFit'      subsurface light value of fitted data
% 'darkValue'            tag-averaged dark PAR value (median of PAR below dark depth)
% 'PAR_QC'               PAR QC flag for the full profile (1: good profile, 2: profile contains any bad values)
% 'SaturationDepth'      depth above which light sensor maxes out
% 'SaturationValue'      PAR value at which light sensor maxes out
% 'quenchDepth'          depth of quenching threshold 15 umol m-2 s-1 (Xing et al. 2018)
% 'Zeu'                  euphotic depth (Morel and Berthon, 1989)
% 'ZeuInterp'            interpolated/smoothed euphotic depths
% 'FirstOD'              First optical depth
% 'meanKd'               mean attenuation coefficient Kd (calculated where LIGHT is non DARK non SAT)
% 'attSlope'             attenuation slope
% 'attSlopeBound'        attenuation slope between pre-defined upper and lower depth bound
% 'predIntChla'          prediction of integrated CHLA from light attenuation slope linear regression
% 
% INPUT ARGUMENTS
% Data - Seal tag/BGC-Argo raw data, processed data and metadata
%   structure, created in loadData.m
% ProfileInfo - Profile specific information/metadata (general, PAR, IRR490, FLUO)
%   structure, created in loadData.m
% defaultPars - default processing parameters
%   structure, created in defaultPars.m
% dataType - Type of radiometric data
%   string, "LIGHT" (for seal tags), "PAR" or "IRR490" (for floats)
%
% OUTPUT
% Data - Seal tag/BGC-Argo raw data, processed data and metadata
%   structure, processed data stored after each processing step
% ProfileInfo - Profile specific information/metadata (general, PAR, IRR490, FLUO)
%   structure, profile info listed in table format

%% CMD message: start
if isempty(getCurrentTask)
    switch dataType
        case 'LIGHT'
            fprintf('Processing <strong>LIGHT</strong> data...');
        case 'PAR'
            fprintf('Processing <strong>PAR</strong> data...');
        case 'IRR490'
            fprintf('Processing <strong>downwelling irradiance at 490nm</strong> data...');
    end
end

%% Create parData info table
var_names = {...
    'Profile', ...              % profile number
    'noData',...                % no data in profile
    'surfaceValue', ...         % subsurface light value (first non NaN value of light)
    'surfaceValueFit', ...      % subsurface light value of fitted data
    'darkValue', ...            % tag-averaged dark PAR value (median of PAR below dark depth)
    'profileQC', ...               % Profile QC flag for the full profile (1: good profile, 2: profile contains any bad values)
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
ProfileInfo.(dataType) = array2table(NaN(Data.MetaData.nProfs, numel(var_names)),'VariableNames',var_names);

%% Populate initial information
% Profile #
ProfileInfo.(dataType).Profile = (1 : Data.MetaData.nProfs)';

% Set profiles with constant PAR values NaN (due to sensor issues)
Data.Processed.(dataType).log.Reg(:,range(Data.Processed.(dataType).log.Reg) == 0) = NaN;
Data.Processed.(dataType).lin.Reg(:,range(Data.Processed.(dataType).lin.Reg) == 0) = NaN;

% Find profiles without usable data (all NaNs)
ProfileInfo.(dataType).noData = all(isnan(Data.Processed.(dataType).log.Reg))' | range(Data.Processed.(dataType).log.Reg)' == 0;

% Shallowest and deepest available observation
firstObs = find_ndim(isfinite(Data.Processed.(dataType).log.Reg),1,'first')';
lastObs = find_ndim(isfinite(Data.Processed.(dataType).log.Reg),1,'last')';

% subsurface light value (first non-NaN value)
nonZeroIndices = find(firstObs(:,1) ~= 0);
ProfileInfo.(dataType).surfaceValue(nonZeroIndices,1) = arrayfun(@(x) Data.Processed.(dataType).lin.Reg(firstObs(x,1), x), nonZeroIndices);

%% Processing: Dark correction, sensor saturation removal, spline fitting
% Initiate processed data fields
Data.Processed.(dataType).log.RegDrk = NaN(size(Data.Processed.DEPTH));
Data.Processed.(dataType).lin.RegDrk = NaN(size(Data.Processed.DEPTH));
Data.Processed.(dataType).log.RegDrkSat = NaN(size(Data.Processed.DEPTH));
Data.Processed.(dataType).lin.RegDrkSat = NaN(size(Data.Processed.DEPTH));
Data.Processed.(dataType).log.RegDrkSatFitAll = NaN(size(Data.Processed.DEPTH));
Data.Processed.(dataType).lin.RegDrkSatFitAll = NaN(size(Data.Processed.DEPTH));
Data.Processed.(dataType).log.RegDrkSatFitSurf = NaN(size(Data.Processed.DEPTH));
Data.Processed.(dataType).lin.RegDrkSatFitSurf = NaN(size(Data.Processed.DEPTH));
Data.Processed.(dataType).log.RegDrkSatFitSurfSmooth = NaN(size(Data.Processed.DEPTH));
Data.Processed.(dataType).lin.RegDrkSatFitSurfSmooth = NaN(size(Data.Processed.DEPTH));
Data.Processed.(dataType).log.RegDrkSatFitBnd = NaN(size(Data.Processed.DEPTH));
Data.Processed.(dataType).lin.RegDrkSatFitBnd = NaN(size(Data.Processed.DEPTH));
Data.Processed.(dataType).QC = NaN(size(Data.Processed.DEPTH));
Data.Processed.(dataType).Kd.FitAll = NaN(size(Data.Processed.(dataType).log.RegDrkSat));
Data.Processed.(dataType).Kd.FitBnd = NaN(size(Data.Processed.DEPTH));
Data.Processed.CHL_LFM = NaN(size(Data.Processed.(dataType).log.RegDrkSat));

%% Only proceed if there is any usable data
if ~all(ProfileInfo.(dataType).noData)
    %% Dark correction
    % Dark values are calculated from a random subsample of deep profiles (>100 m, N = 10 or as defined in
    % defaultPars.PAR.nDarkProfiles) and used to calculate a tag-specific dark value. The randomised subsampling of deep
    % profiles is repeated nRep times to yield a more robust estimate. The tag-specific dark value is calculated as the
    % median of the median of all individual dark values (per subsample). NB: The random number generator algorithm is reset
    % within the script to ensure that the results do not change when reprocessing data.
    % The calculation of dark values from PAR profiles is based on a Lilliefors normality test (following Organelli et al.,
    % 2016). The test is iteratively applied to PAR values from the deepest available observation (z0) to gradually shallower
    % depths (z1). While the null hypothesis is rejected (normal distribution of values between z0 and z1), PAR is assumed
    % to be part of the dark signal. When the null-hypothesis can no longer be rejected, z1 is taken to be the upper end of
    % the dark signal. The median of all PAR values between z0 and z1 is the profile specific dark value.

    PAR_lin_RegDrk = Data.Processed.(dataType).lin.Reg;
    PAR_log_RegDrk = Data.Processed.(dataType).log.Reg;

    rng(0,'twister')                                % Reset random number generator algorithm (for reproducibility of the random subsampling)
    deepProfs = find(lastObs > 100)';               % Indices of all deep profiles (deeper than 100 m)
    nRep = 10;                                      % Number of repetitions of random deep profile subsamples
    darkValues = NaN(Data.MetaData.nProfs,nRep);    % Dark value matrix
    darkDepths = NaN(Data.MetaData.nProfs,nRep);    % Dark depth matrix
    for iRep = 1 : nRep
        % Randomly subsample deep profiles to be used for the dark value calculation (N defined in setDefaults)
        nSamples = min(defaultPars.PAR.nDarkProfiles,numel(deepProfs)); % N subsamples cannot be less than the number of deep profiles to sample from
        % iDarkProfiles = datasample(deepProfs,nSamples,'Replace',false); % Indices of subsampled profiles
        iDarkProfiles = datasample(deepProfs,nSamples,'Replace',false); % Indices of subsampled profiles

        % Loop through profiles selected for dark value calculation
        for iP = iDarkProfiles
            %% 1) Detect and remove constant part of the PAR profile to avoid normality test on constant data series
            % See this example:
            % [h_temp,p_temp] = lillietest(ones(1,defVals.depthInterpGrid(end)),'Alpha',0.01,'Distribution','normal');
            % returns h_temp = 1, i.e. null hypothesis not rejected (= distribution is normal)

            % iteratively move down the water column and check the range of PAR values between depth iZ and the deepest PAR value

            for iZ = firstObs(iP) : lastObs(iP) - 1

                iZ_PAR_range = range(PAR_lin_RegDrk(iZ:lastObs(iP),iP));

                if iZ_PAR_range <= defaultPars.PAR.delta_constantPAR
                    Z_const = iZ + 1;           % depth of first value belonging to constant signal
                    PAR_lin_RegDrk(Z_const:end,iP) = NaN;
                    PAR_log_RegDrk(Z_const:end,iP) = NaN;
                    lastObs(iP) = Z_const-1;    % Update last available observation depth

                    break   % end for-loop when PAR range falls below threshold value where PAR is considered constant
                end

            end

            %% 2) Detect PAR dark value (Organelli et al., 2016)
            % Initialize values
            PAR_dark = NaN;     % dark PAR value
            Z_dark = NaN;       % dark depth

            % Lilliefors test requires at least 4 valid observations
            if lastObs(iP) - firstObs(iP) + 1 >= 4
                % counter and required number of consecutive accepted null hypothesis tests determining the end of the dark signal
                ct_h1 = 0;
                n_h1 = 3;
                % iteratively move up the water column and check normality of distribution of PAR values from the deepest value to depth i
                for iZ = lastObs(iP) - 3 : -1 : firstObs(iP)
                    % Warning: P is less than the smallest tabulated value, returning 0.001.
                    warning('off','stats:lillietest:OutOfRangePLow')
                    warning('off','stats:lillietest:OutOfRangePHigh')

                    % Lilliefor normality test
                    [h_lillie,~] = lillietest(PAR_lin_RegDrk(iZ:lastObs(iP),iP),'Alpha',0.01,'Distribution','normal');

                    if h_lillie == 1
                        % if null hypothesis is not rejected (non-normal distribution above dark signal) increase counter by 1
                        ct_h1 = ct_h1 + 1;

                        if ct_h1 == n_h1 || iZ == firstObs(iP)
                            % If non-dark signal is found at n consecutive depths (or if dark signal extends to the surface)
                            % >> dark signal starts at iZ+n
                            Z_dark = iZ + n_h1;
                            PAR_dark = median(PAR_lin_RegDrk(Z_dark:lastObs(iP),iP));

                            break % end for-loop
                        end
                    else
                        ct_h1 = 0;
                    end
                end
            end

            % Write to parData table
            darkValues(iP,iRep) = PAR_dark;	 % all dark PAR values of the selected profiles
            % Write to parData table
            darkDepths(iP,iRep) = Z_dark;	 % all dark PAR values of the selected profiles

        end
    end

    % Calculate median of calulated dark values and apply to every profile
    ProfileInfo.(dataType).darkValue(:) = median(median(darkValues,1,'omitnan'),2,'omitnan');

    % Offset all non-dark PAR values by the dark value
    PAR_log_RegDrk = PAR_log_RegDrk - log(ProfileInfo.(dataType).darkValue(1));

    % Remove negative/0 PAR values and linearly interpolate gaps
    PAR_log_RegDrk(PAR_lin_RegDrk<=0) = NaN;
    PAR_log_RegDrk = fillmissing(PAR_log_RegDrk,'linear',1,'EndValues','none');

    % Write to platform_processed structure
    Data.Processed.(dataType).log.RegDrk = PAR_log_RegDrk;
    Data.Processed.(dataType).lin.RegDrk = exp(PAR_log_RegDrk);

    % Update last non-NaN parData table column
    firstObs = find_ndim(isfinite(Data.Processed.(dataType).log.RegDrk),1,'first')';
    lastObs = find_ndim(isfinite(Data.Processed.(dataType).log.RegDrk),1,'last')';

    %% Attenuation slope (based on logarithmic PAR values)
    % Slope calculated over non-dark section of PAR profile
    for iP = 1 : Data.MetaData.nProfs
        if firstObs(iP) > 0
            atten_slope = ...
                (Data.Processed.(dataType).log.RegDrk(lastObs(iP),iP) - Data.Processed.(dataType).log.RegDrk(firstObs(iP),iP)) /...
                (defaultPars.depthInterpGrid(lastObs(iP)) - defaultPars.depthInterpGrid(firstObs(iP)));

            % write parData table
            ProfileInfo.(dataType).attSlope(iP)        = atten_slope;	% slope of attenuation (calculated where PAR is non DARK)
        end
    end

    %% QC surface light values
    % Assign QC flags to PAR profiles and individual observations based on how closely light follows the expected exponential
    % increasing trend towards the surface. Sensor saturation is assumed where consecutive near-surface observations reach a
    % plateau that is within the highest values measured by the respective tag.
    % NB: Saturation correction is only applied to seal tag data (BGC-Argo PAR sensors do not saturate)
    %
    % The following QC flags are assigned:
    % 0 - Profile/observations that do not exceed 1 umol quanta m-2 s-1 (no QC performed)
    % 1 - Profile/observations where light increases exponentially towards the surface without significant departures
    % 2 - Profile/observations where light diverges significantly from the expected exponential increase towards the surface
    % 3 - Profile/observations where unexpected changes in light at the surface are suspected to be due to light sensor saturation
    % NaN - Profile without light data (following dark correction)

    ProfileInfo.(dataType).profileQC = NaN(Data.MetaData.nProfs,1);
    Data.Processed.(dataType).log.RegDrkSat = Data.Processed.(dataType).log.RegDrk;
    saturationPoints = NaN(Data.MetaData.nProfs,1);

    % Define a threshold PAR value above which observations will be considered for QC
    highPAR_thresh = 0; % log(PAR) = 0 -> PAR = 1

    % Near-maximum PAR value across all profiles (defined as percentile), considering only PAR values above the defined
    % threshold to not be biased by negligible light values at depth. Used as a reference to determine if constant light
    % signal at the surface is saturated or cloud/ice-shading related.
    tagHighPAR = prctile(Data.Processed.(dataType).log.RegDrk(Data.Processed.(dataType).log.RegDrk>=highPAR_thresh),95,[1,2]);

    for iP = 1 : Data.MetaData.nProfs
        % Isolate PAR profile data for QC
        logPAR = Data.Processed.(dataType).log.RegDrkSat(:,iP);
        
        if any(isfinite(logPAR))
            % Set QC flags to 0/1 where PAR is above/below the defined threshold
            Data.Processed.(dataType).QC(logPAR<highPAR_thresh,iP) = 0;
            Data.Processed.(dataType).QC(logPAR>=highPAR_thresh,iP) = 1;

            if any(logPAR>=highPAR_thresh)
                % If any PAR observations are above the threshold, assign profile QC flag 1
                ProfileInfo.(dataType).profileQC(iP,1) = 1;
                
                % Remove observations where PAR is below the threshold (not considered in further QC)
                logPAR(logPAR<0) = NaN;

                % Mean and SD of the derivative of PAR above the defined threshold (i.e. light attenuation)
                derivPAR_mean = mean(diff(logPAR),'omitnan');
                derivPAR_std = std(diff(logPAR),0,'omitnan');

                % Lower threshold of expected light attenuation between depth levels
                % >> mean+SD (must not be positive, i.e. light has to decrease with depth)
                derivPAR_thresh = min([derivPAR_mean+derivPAR_std,0]);

                % Flag observations where the derivative of PAR exceeds the threshold value, i.e. PAR values that do not
                % increase towards the surface as expected (either taper out or decrease)
                flagged = find([diff(logPAR)>=derivPAR_thresh ;false]);

                if ~isempty(flagged)
                    % Assign QC flag 2 to flagged observations and the profile
                    Data.Processed.(dataType).QC(flagged,iP) = 2;
                    ProfileInfo.(dataType).profileQC(iP,1) = 2;

                    % Proceed to check for potential saturation signal:
                    % Find consecutive flagged observations near the surface. Consecutive flagged observations can be
                    % interrupted by up to 5 unflagged observations to account for unstable/variable saturation signal.
                    flagged_shallow_obs = flagged(1) : flagged(find(diff([flagged;inf]) > 5,1,'first'));

                    % These observations are assumed to be saturated if ...
                    % (1) they start at or near the surface (within the 3 shallowest available observations)
                    % (2) they cover at least 3 observations
                    % (3) their mean value is within >98% of the maximum light value measured across all profiles
                    % (4) their differential is less than or equal to 1 (no decline in PAR)
                    % (5) the profile has not been flagged as being potentially sea ice covered
                    %     --> NB: flagged obs are by definition (near-)constant
                    isSaturated = ...
                        flagged_shallow_obs(1) <= firstObs(iP)+3 & ...
                        numel(flagged_shallow_obs) >= 3 & ...
                        median(logPAR(flagged_shallow_obs)) >= tagHighPAR & ...
                        median(diff(logPAR(flagged_shallow_obs))) - mad(diff(logPAR(flagged_shallow_obs)),1) <= 0 & ...
                        ~ProfileInfo.General.SeaIceCovered(iP);

                    % Saturation correction is only applied to seal tag data (BGC-Argo PAR sensors do not saturate)
                    if isSaturated && strcmp(Data.MetaData.platform_type,'sealtag')
                        % Extend depths identified as being saturated to the first available observation
                        flagged_shallow_obs = firstObs(iP) : flagged_shallow_obs(end);
                        
                        % Assign QC flag 3 to flagged shallow observations and profile
                        Data.Processed.(dataType).QC(flagged_shallow_obs,iP) = 3;
                        saturationPoints(iP) = mean(logPAR(flagged_shallow_obs));
                        ProfileInfo.(dataType).profileQC(iP,1) = 3;

                        % Remove saturated values
                        Data.Processed.(dataType).log.RegDrkSat(flagged_shallow_obs,iP) = NaN;

                    end
                end
            else
                % If all PAR observations are below the threshold, assign profile QC flag 0 (no further QC performed)
                ProfileInfo.(dataType).profileQC(iP,1) = 0;
            end
        end
    end

    % Save saturation points to profile PAR info table
    ProfileInfo.(dataType).SaturationValue = exp(saturationPoints);
    
    % Write linear PAR data to platform_processed structure
    Data.Processed.(dataType).lin.RegDrkSat = exp(Data.Processed.(dataType).log.RegDrkSat);

    %% Spline/log fit to observations
    Data = processRadiometry_fit(Data,ProfileInfo,defaultPars,dataType);

    %% Final things
    % - New subsurface PAR value based on fitted data
    % - Euphotic depth (1% of the surface PAR depth, after Morel and Berthon, 1989)
    iSurfVal_fit = find_ndim(isfinite(Data.Processed.(dataType).lin.RegDrkSatFitAll),1,'first');
    finiteIndices = find(iSurfVal_fit>0);
    for iP = finiteIndices
        ProfileInfo.(dataType).surfaceValueFit(iP,1) = Data.Processed.(dataType).lin.RegDrkSatFitAll(iSurfVal_fit(iP),iP);
        Zeu = -find(Data.Processed.(dataType).lin.RegDrkSatFitAll(:,iP) / ProfileInfo.(dataType).surfaceValueFit(iP) < 0.01,1,'first');
        if ~isempty(Zeu)
            ProfileInfo.(dataType).Zeu(iP,1) = Zeu;
        end
    end
    % fill missing and smooth Zeu data (window size: ~n profiles per day)
    ProfileInfo.(dataType).ZeuInterp = fillmissing(ProfileInfo.(dataType).Zeu,'linear','EndValues','none');
    ProfileInfo.(dataType).ZeuInterp = movmedian(ProfileInfo.(dataType).ZeuInterp,ceil(Data.MetaData.nProfs / max(ProfileInfo.General.DeployDay)),'omitnan','Endpoints','shrink');

    % penetration depth (first optical depth)
    ProfileInfo.(dataType).FirstOD = ProfileInfo.(dataType).ZeuInterp / 4.6;

    % depth of quenching threshold 15 umol m-2 s-1 (Xing et al. 2018)
    iZ_quench = find_ndim(Data.Processed.(dataType).lin.RegDrkSatFitAll > 15, 1, 'last')';
    iZ_quench(iZ_quench == 0) = NaN;
    ProfileInfo.(dataType).quenchDepth = -iZ_quench;

    % (ML-)mean Kd
    ProfileInfo.(dataType).meanKd = mean(Data.Processed.(dataType).Kd.FitAll,1,'omitnan')';
    for iP = 1 : Data.MetaData.nProfs
        iZ_MLD = abs(round(ProfileInfo.General.MLD(iP)));
        if isfinite(iZ_MLD)
            ProfileInfo.(dataType).meanKdML(iP) = mean(Data.Processed.(dataType).Kd.FitAll(1:iZ_MLD,iP),1,'omitnan');
        end
    end
end

%% CMD message: done
if isempty(getCurrentTask)
    fprintf('\b\b \x2713\n')
end

end