function [tagProcessed,ProfileInfo_PAR] = processPAR(tagMetadata,tagProcessed,ProfileInfo,defaultPars)
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
    'darkValue', ...            % tag-averaged dark PAR value (median of PAR below dark depth)
    'PAR_QC', ...               % PAR QC flag for the full profile (1: good profile, 2: profile contains any bad values)
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
ProfileInfo_PAR = array2table(NaN(tagMetadata.nProfs, numel(var_names)),'VariableNames',var_names);

%% Populate initial information
% Profile #
ProfileInfo_PAR.Profile = (1 : tagMetadata.nProfs)';

% Set profiles with constant PAR values NaN (due to sensor issues)
tagProcessed.PAR.log.Reg(:,range(tagProcessed.PAR.log.Reg) == 0) = NaN;
tagProcessed.PAR.lin.Reg(:,range(tagProcessed.PAR.lin.Reg) == 0) = NaN;

% Find profiles without usable data (all NaNs)
ProfileInfo_PAR.noData = all(isnan(tagProcessed.PAR.log.Reg))' | range(tagProcessed.PAR.log.Reg)' == 0;

% Shallowest and deepest available observation
firstObs = find_ndim(isfinite(tagProcessed.PAR.log.Reg),1,'first')';
lastObs = find_ndim(isfinite(tagProcessed.PAR.log.Reg),1,'last')';

% subsurface light value (first non-NaN value)
nonZeroIndices = find(firstObs(:,1) ~= 0);
ProfileInfo_PAR.surfaceValue(nonZeroIndices,1) = arrayfun(@(x) tagProcessed.PAR.lin.Reg(firstObs(x,1), x), nonZeroIndices);

%% Processing: Dark correction, sensor saturation removal, spline fitting
% Initiate processed data fields
tagProcessed.PAR.log.RegDrk = NaN(size(tagProcessed.DEPTH));
tagProcessed.PAR.lin.RegDrk = NaN(size(tagProcessed.DEPTH));
tagProcessed.PAR.log.RegDrkSat = NaN(size(tagProcessed.DEPTH));
tagProcessed.PAR.lin.RegDrkSat = NaN(size(tagProcessed.DEPTH));
tagProcessed.PAR.log.RegDrkSatFitAll = NaN(size(tagProcessed.DEPTH));
tagProcessed.PAR.lin.RegDrkSatFitAll = NaN(size(tagProcessed.DEPTH));
tagProcessed.PAR.log.RegDrkSatFitSurf = NaN(size(tagProcessed.DEPTH));
tagProcessed.PAR.lin.RegDrkSatFitSurf = NaN(size(tagProcessed.DEPTH));
tagProcessed.PAR.log.RegDrkSatFitSurfSmooth = NaN(size(tagProcessed.DEPTH));
tagProcessed.PAR.lin.RegDrkSatFitSurfSmooth = NaN(size(tagProcessed.DEPTH));
tagProcessed.PAR.log.RegDrkSatFitBnd = NaN(size(tagProcessed.DEPTH));
tagProcessed.PAR.lin.RegDrkSatFitBnd = NaN(size(tagProcessed.DEPTH));
tagProcessed.PAR.QC = NaN(size(tagProcessed.DEPTH));
tagProcessed.PAR.Kd.FitAll = NaN(size(tagProcessed.PAR.log.RegDrkSat));
tagProcessed.PAR.Kd.FitBnd = NaN(size(tagProcessed.DEPTH));
tagProcessed.CHL_LFM = NaN(size(tagProcessed.PAR.log.RegDrkSat));

%% Only proceed if there is any usable data
if ~all(ProfileInfo_PAR.noData)
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

    PAR_lin_RegDrk = tagProcessed.PAR.lin.Reg;
    PAR_log_RegDrk = tagProcessed.PAR.log.Reg;

    rng(0,'twister')                            % Reset random number generator algorithm (for reproducibility of the random subsampling)
    deepProfs = find(lastObs > 100)';           % Indicies of all deep profiles (deeper than 100 m)
    nRep = 10;                                  % Number of repetitions of random deep profile subsamples
    darkValues = NaN(tagMetadata.nProfs,nRep);  % Dark value matrix
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
            % returns h_temp = 1, i.e. rejection of the null hypothesis (= distribution is normal)

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

                        if ct_h1 == n_h1
                            % If non-dark signal is found at n consecutive depths then the dark signal starts at iZ+n
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
        end
    end

    % Calculate median of calulated dark values and apply to every profile
    ProfileInfo_PAR.darkValue(:) = median(median(darkValues,1,'omitnan'),2,'omitnan');

    % Offset all non-dark PAR values by the dark value
    PAR_lin_RegDrk = PAR_lin_RegDrk - ProfileInfo_PAR.darkValue(1);

    % Remove negative/0 PAR values and linearly interpolate gaps
    PAR_lin_RegDrk(PAR_lin_RegDrk<=0) = NaN;
    PAR_lin_RegDrk = fillmissing(PAR_lin_RegDrk,'linear',1,'EndValues','none');

    % Write to platform_processed structure
    tagProcessed.PAR.lin.RegDrk = PAR_lin_RegDrk;
    tagProcessed.PAR.log.RegDrk = log(PAR_lin_RegDrk);

    % Update last non-NaN parData table column
    firstObs = find_ndim(isfinite(tagProcessed.PAR.log.RegDrk),1,'first')';
    lastObs = find_ndim(isfinite(tagProcessed.PAR.log.RegDrk),1,'last')';

    %% Attenuation slope (based on logarithmic PAR values)
    % Slope calculated over non-dark section of PAR profile
    for iP = 1 : tagMetadata.nProfs
        if firstObs(iP) > 0
            atten_slope = ...
                (tagProcessed.PAR.log.RegDrk(lastObs(iP),iP) - tagProcessed.PAR.log.RegDrk(firstObs(iP),iP)) /...
                (defaultPars.depthInterpGrid(lastObs(iP)) - defaultPars.depthInterpGrid(firstObs(iP)));

            % write parData table
            ProfileInfo_PAR.attSlope(iP)        = atten_slope;	% slope of attenuation (calculated where PAR is non DARK)
        end
    end

    %% QC surface light values
    % Assign QC flags to PAR profiles and individual observations based on how closely light follows the expected exponential
    % increasing trend towards the surface. Sensor saturation is assumed where consecutive near-surface observations reach a
    % plateau that is within the highest values measured by the respective tag.
    %
    % The following QC flags are assigned:
    % 0 - Profile/observations that do not exceed 1 umol quanta m-2 s-1 (no QC performed)
    % 1 - Profile/observations where light increases exponentially towards the surface without significant departures
    % 2 - Profile/observations where light diverges significantly from the expected exponential increase towards the surface
    % 3 - Profile/observations where unexpected changes in light at the surface are suspected to be due to light sensor saturation
    % NaN - Profile without light data (following dark correction)

    ProfileInfo_PAR.PAR_QC = NaN(tagMetadata.nProfs,1);
    tagProcessed.PAR.log.RegDrkSat = tagProcessed.PAR.log.RegDrk;
    saturationPoints = NaN(tagMetadata.nProfs,1);

    % Define a threshold PAR value above which observations will be considered for QC
    highPAR_thresh = 0; % log(PAR) = 0 -> PAR = 1

    % Near-maximum PAR value across all profiles (defined as percentile), considering only PAR values above the defined
    % threshold to not be biased by negligible light values at depth. Used as a reference to determine if constant light
    % signal at the surface is saturated or cloud/ice-shading related.
    tagHighPAR = prctile(tagProcessed.PAR.log.RegDrk(tagProcessed.PAR.log.RegDrk>=highPAR_thresh),95,[1,2]);

    for iP = 1 : tagMetadata.nProfs
        % Isolate PAR profile data for QC
        logPAR = tagProcessed.PAR.log.RegDrkSat(:,iP);
        
        if any(isfinite(logPAR))
            % Set QC flags to 0/1 where PAR is above/below the defined threshold
            tagProcessed.PAR.QC(logPAR<highPAR_thresh,iP) = 0;
            tagProcessed.PAR.QC(logPAR>=highPAR_thresh,iP) = 1;

            if any(logPAR>=highPAR_thresh)
                % If any PAR observations are above the threshold, assign profile QC flag 1
                ProfileInfo_PAR.PAR_QC(iP,1) = 1;
                
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
                    tagProcessed.PAR.QC(flagged,iP) = 2;
                    ProfileInfo_PAR.PAR_QC(iP,1) = 2;

                    % Proceed to check for potential saturation signal:
                    % Find consecutive flagged observations near the surface. Consecutive flagged observations can be
                    % interrupted by up to 5 unflagged observations to account for unstable/variable saturation signal.
                    flagged_shallow_obs = flagged(1) : flagged(find(diff([flagged;inf]) > 5,1,'first'));

                    % These observations are assumed to be saturated if ...
                    % (1) they start at or near the surface (within the 3 shallowest available observations)
                    % (2) they cover at least 3 observations
                    % (3) their mean value is within >98% of the maximum light value measured across all profiles
                    % (4) their differential is less than or equal to 1 (no decline in PAR)
                    %     --> NB: flagged obs are by definition (near-)constant
                    isSaturated = ...
                        flagged_shallow_obs(1) <= firstObs(iP)+3 & ...
                        numel(flagged_shallow_obs) >= 3 & ...
                        median(logPAR(flagged_shallow_obs)) >= tagHighPAR & ...
                        median(diff(logPAR(flagged_shallow_obs))) - mad(diff(logPAR(flagged_shallow_obs)),1) <= 0;

                    if isSaturated
                        % Extend depths identified as being saturated to the first available observation
                        flagged_shallow_obs = firstObs(iP) : flagged_shallow_obs(end);
                        
                        % Assign QC flag 3 to flagged shallow observations and profile
                        tagProcessed.PAR.QC(flagged_shallow_obs,iP) = 3;
                        saturationPoints(iP) = mean(logPAR(flagged_shallow_obs));
                        ProfileInfo_PAR.PAR_QC(iP,1) = 3;

                        % Remove saturated values
                        tagProcessed.PAR.log.RegDrkSat(flagged_shallow_obs,iP) = NaN;

                    end
                end
            else
                % If all PAR observations are below the threshold, assign profile QC flag 0 (no further QC performed)
                ProfileInfo_PAR.PAR_QC(iP,1) = 0;
            end
        end
    end

    % Save saturation points to profile PAR info table
    ProfileInfo_PAR.SaturationValue = exp(saturationPoints);
    
    % Write linear PAR data to platform_processed structure
    tagProcessed.PAR.lin.RegDrkSat = exp(tagProcessed.PAR.log.RegDrkSat);

    %% Spline/log fit to observations
    tagProcessed = processPAR_fit(tagProcessed,ProfileInfo,ProfileInfo_PAR,defaultPars);

    %% Final things
    % - New subsurface PAR value based on fitted data
    % - Euphotic depth (1% of the surface PAR depth, after Morel and Berthon, 1989)
    iSurfVal_fit = find_ndim(isfinite(tagProcessed.PAR.lin.RegDrkSatFitAll),1,'first');
    finiteIndices = find(iSurfVal_fit>0);
    for iP = finiteIndices
        ProfileInfo_PAR.surfaceValueFit(iP,1) = tagProcessed.PAR.lin.RegDrkSatFitAll(iSurfVal_fit(iP),iP);
        Zeu = -find(tagProcessed.PAR.lin.RegDrkSatFitAll(:,iP) / ProfileInfo_PAR.surfaceValueFit(iP) < 0.01,1,'first');
        if ~isempty(Zeu)
            ProfileInfo_PAR.Zeu(iP,1) = Zeu;
        end
    end
    % fill missing and smooth Zeu data (window size: ~n profiles per day)
    ProfileInfo_PAR.ZeuInterp = fillmissing(ProfileInfo_PAR.Zeu,'linear','EndValues','none');
    ProfileInfo_PAR.ZeuInterp = movmedian(ProfileInfo_PAR.ZeuInterp,ceil(tagMetadata.nProfs / max(ProfileInfo.DeployDay)),'omitnan','Endpoints','shrink');

    % penetration depth (first optical depth)
    ProfileInfo_PAR.FirstOD = ProfileInfo_PAR.ZeuInterp / 4.6;

    % depth of quenching threshold 15 umol m-2 s-1 (Xing et al. 2018)
    iZ_quench = find_ndim(tagProcessed.PAR.lin.RegDrkSatFitAll > 15, 1, 'last')';
    iZ_quench(iZ_quench == 0) = NaN;
    ProfileInfo_PAR.quenchDepth = -iZ_quench;

    % (ML-)mean Kd
    ProfileInfo_PAR.meanKd = mean(tagProcessed.PAR.Kd.FitAll,1,'omitnan')';
    for iP = 1 : tagMetadata.nProfs
        iZ_MLD = abs(round(ProfileInfo.MLD(iP)));
        if isfinite(iZ_MLD)
            ProfileInfo_PAR.meanKdML(iP) = mean(tagProcessed.PAR.Kd.FitAll(1:iZ_MLD,iP),1,'omitnan');
        end
    end
end

%% CMD message: done
fprintf('\b\b \x2713\n')

end