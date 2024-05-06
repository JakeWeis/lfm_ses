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
tagProcessed.PAR.log_Reg(:,range(tagProcessed.PAR.log_Reg) == 0) = NaN;
tagProcessed.PAR.lin_Reg(:,range(tagProcessed.PAR.lin_Reg) == 0) = NaN;

% Find profiles without usable data (all NaNs)
ProfileInfo_PAR.noData = all(isnan(tagProcessed.PAR.log_Reg))' | range(tagProcessed.PAR.log_Reg)' == 0;

% Shallowest and deepest available observation
firstObs = find_ndim(isfinite(tagProcessed.PAR.log_Reg),1,'first')';
lastObs = find_ndim(isfinite(tagProcessed.PAR.log_Reg),1,'last')';

% subsurface light value (first non-NaN value)
nonZeroIndices = find(firstObs(:,1) ~= 0);
ProfileInfo_PAR.surfaceValue(nonZeroIndices,1) = arrayfun(@(x) tagProcessed.PAR.lin_Reg(firstObs(x,1), x), nonZeroIndices);

% Only proceed with processing if there is any usable data
if ~all(ProfileInfo_PAR.noData)
    %% Dark depth/value & attenuation slope calculation
    PAR_lin_RegDrk = tagProcessed.PAR.lin_Reg;
    PAR_log_RegDrk = tagProcessed.PAR.log_Reg;

    % Calculate dark value from the n deepest profiles >100 m (n = defaultPars.PAR.nDarkProfiles)
    [nMaxDepths,iMaxDepthProfs] = maxk(lastObs,defaultPars.PAR.nDarkProfiles);
    iDarkProfiles = iMaxDepthProfs(nMaxDepths > 100)';

    darkDepth = NaN(tagMetadata.nProfs,1);

    % Loop through profiles selected for processing
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
        darkDepth(iP) = -Z_dark;     % first dark depth (index, depth)
        ProfileInfo_PAR.darkValue(iP) = PAR_dark;	% dark PAR value

        %% 3) Attenuation slope (based on logarithmic PAR values)
        % Slope calculated over non-dark section of PAR profile
        if isfinite(Z_dark) && firstObs(iP) > 0
            atten_slope = ...
                (PAR_log_RegDrk(Z_dark-1,iP) - PAR_log_RegDrk(firstObs(iP),iP)) /...
                (defaultPars.depthInterpGrid(Z_dark-1) - defaultPars.depthInterpGrid(firstObs(iP)));

            % Slope calculated between defined upper and lower depth bound
            % (only if dark depth is deeper than lower bound)
            if defaultPars.lowerDepthBound > darkDepth(iP)
                atten_slope_PART = ...
                    (PAR_log_RegDrk(abs(defaultPars.lowerDepthBound),iP) - PAR_log_RegDrk(firstObs(iP),iP)) /...
                    (tagProcessed.DEPTH(abs(defaultPars.lowerDepthBound),iP) - tagProcessed.DEPTH(firstObs(iP),iP));
            else
                atten_slope_PART = NaN;
            end

            % write parData table
            ProfileInfo_PAR.attSlope(iP)        = atten_slope;	% slope of attenuation (calculated where PAR is non DARK)
            ProfileInfo_PAR.attSlopeBound(iP)    = atten_slope_PART; % slope of attenuation calculated on [upperDepthBound:lowerDepthBound] interval
        end
    end

    %% Dark correction
    % Calculate median of calulated dark values and apply to every profile
    ProfileInfo_PAR.darkValue(:) = median(ProfileInfo_PAR.darkValue,'omitnan');

    % Offset all non-dark PAR values by the dark value
    PAR_lin_RegDrk = PAR_lin_RegDrk - ProfileInfo_PAR.darkValue(1);

    % Remove negative/0 PAR values and linearly interpolate gaps
    PAR_lin_RegDrk(PAR_lin_RegDrk<=0) = NaN;
    PAR_lin_RegDrk = fillmissing(PAR_lin_RegDrk,'linear',1,'EndValues','none');

    % Write to platform_processed structure
    tagProcessed.PAR.lin_RegDrk = PAR_lin_RegDrk;
    tagProcessed.PAR.log_RegDrk = log(PAR_lin_RegDrk);

    % Update last non-NaN parData table column
    lastObs = find_ndim(isfinite(tagProcessed.PAR.log_RegDrk),1,'last')';

    %% Saturation depth
    % PAR_lin_RegDrkSat = tagProcessed.PAR.lin_RegDrk;
    % 
    % for iP = 1 : tagMetadata.nProfs
    %     % derivative of log dark-corrected PAR
    %     diffPARlog = diff(PAR_log_RegDrk(:,iP));
    %     % mean of derivative
    %     meandiffPARlog = mean(diffPARlog,'omitnan');
    % 
    %     % check saturation depth only if profile starts with descrease>mean derivative
    %     % else no saturation depth and saturation value is defined to maximum PAR value measured
    %     if isfinite(meandiffPARlog) && diffPARlog(firstObs(iP)) > meandiffPARlog
    %         Z_saturation = firstObs(iP) + 1;
    %         while diffPARlog(Z_saturation) >= meandiffPARlog
    %             Z_saturation = Z_saturation + 1;
    %         end
    %         % Remove saturated values
    %         PAR_lin_RegDrkSat(1 : Z_saturation-1,iP) = NaN;
    % 
    %         % Write saturation depth and PAR value to parData table
    %         ProfileInfo_PAR.SaturationDepth(iP) = -Z_saturation;
    %         ProfileInfo_PAR.SaturationValue(iP) = PAR_lin_RegDrkSat(Z_saturation,iP);
    %     end
    % end
    % 
    % % Write to platform_processed structure
    % tagProcessed.PAR.lin_RegDrkSat = PAR_lin_RegDrkSat;
    % tagProcessed.PAR.log_RegDrkSat = log(PAR_lin_RegDrkSat);
    % 
    % % Update first non-NaN parData table column
    % firstObs = find_ndim(isfinite(tagProcessed.PAR.lin_RegDrkSat),1,'first')';

    %% Saturation depth (ATLERNATIVE ID METHOD)
    PAR_lin_RegDrkSat = tagProcessed.PAR.lin_RegDrk;

    % Maximum PAR value in each profile and depths at which maxima occur
    maxPAR = max(PAR_lin_RegDrkSat)';
    i_maxPAR = arrayfun(@(iP) find(PAR_lin_RegDrkSat(:,iP) == maxPAR(iP)), 1 : tagMetadata.nProfs, 'UniformOutput',false)';

    % Identify profiles with apparent saturation points:
    % --> Profiles with maximum PAR values that are constant over 2 or more consecutive depths
    maybeSaturated = find(arrayfun(@(iP) numel(i_maxPAR{iP}) > 1 && ismember(1,diff(i_maxPAR{iP})), 1 : tagMetadata.nProfs));

    % Check for "unsaturated" maxima exceeding the highest apparent saturation maximum:
    % If (a) at least 3 PAR maxima exist that are NOT identified as saturation points and (b) exceed the maximum apparent
    % saturation point on average by >10%.
    % --> This condition identifies cases where saturation is associated with relatively low PAR values and therefore cannot
    % reasonably be assumed to reflect a saturation point.
    lowSaturationPoint = (numel(maxPAR(max(maxPAR(maybeSaturated)) < maxPAR)) > 2 && ...            % (a)
        mean(maxPAR(max(maxPAR(maybeSaturated)) < maxPAR) ./ max(maxPAR(maybeSaturated)))>1.1);     % (b)

    if lowSaturationPoint
        % If the above condition is true, all saturation points are assumed to be apparent only and are therefore ignored.
        isSaturated = [];
    else
        % Else proceed to identify a reasonable range of true saturation values via an iterative approach, determining the
        % largest set of maximum saturation values for which the 20th percentile is within 80% of the 80th percentile (i.e.
        % find a concise set of the highest saturation values). Low saturation values are discarded. This approach assumes
        % that sensor saturation should be more or less consistent across all profiles and can only be reflected by the
        % highest values.
        % Note: Where the initial set of apparent saturation values is already within the specified range, no saturation
        % values are discarded. Else low saturation points are discounted. 

        [~,iSort] = sort(maxPAR(maybeSaturated));   % Apparent saturation values (sorted)
        maybeSaturated = maybeSaturated(iSort);     % Sort apparent saturation indices
        isSaturated = maybeSaturated;
        ct = 0;
        while prctile(maxPAR(isSaturated),20)/prctile(maxPAR(isSaturated),80) < 0.8
            ct = ct + 1;
            isSaturated = maybeSaturated(ct:end);
        end
    end

    for iP = isSaturated
        % Set values above the lowest saturated value NaN
        PAR_lin_RegDrkSat(1:i_maxPAR{iP}(end),iP) = NaN;
        
        % Write saturation depth and PAR value to parData table
        ProfileInfo_PAR.SaturationDepth(iP) = -i_maxPAR{iP}(end);
        ProfileInfo_PAR.SaturationValue(iP) = maxPAR(iP);
    end

    % Write to platform_processed structure
    tagProcessed.PAR.lin_RegDrkSat = PAR_lin_RegDrkSat;
    tagProcessed.PAR.log_RegDrkSat = log(PAR_lin_RegDrkSat);

    % Update first non-NaN parData table column
    firstObs = find_ndim(isfinite(tagProcessed.PAR.lin_RegDrkSat),1,'first')';

    %% Functional fit (full profile)
    tagProcessed.PAR.log_RegDrkSatFitAll = NaN(size(tagProcessed.PAR.log_RegDrkSat));
    tagProcessed.PAR.Kd_all = NaN(size(tagProcessed.PAR.log_RegDrkSat));
    % platform_processed.CHL_LFM.Full = NaN(size(platform_processed.PAR.log_RegDrkSat));

    % Select profiles for fit computation
    lumToFit_all = tagProcessed.PAR.log_RegDrkSat;
    i_profiles = find(...
        ProfileInfo.Processed &...              % Passed processing checks in loadData
        ProfileInfo.Daytime &...                % Daytime profile
        sum(~isnan(tagProcessed.PAR.log_RegDrkSat))' >= 3)';  % At least three values

    % Fit parameters
    % Initialise bspline fit coefficients
    fdCoeffs = zeros(defaultPars.LFM.nBasis,1);
    % Set up functional parameter object
    lfdObj = 1; % original value in script = 2; penalize curvature of acceleration

    % Compute fit
    for iP = i_profiles
        % PAR data to compute fit over
        PAR_fit = lumToFit_all(:,iP);
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
        tagProcessed.PAR.log_RegDrkSatFitAll(fitInterval,iP) = eval_fd(fdFit_PAR,fitInterval);

        % Calculate Kd as the derivative function of the PAR fit
        fdKd_iP = deriv_fd(fdFit_PAR);
        tagProcessed.PAR.Kd_all(fitInterval,iP) = -eval_fd(fdKd_iP,fitInterval);

        %% LFM CHL PREDICTION FROM KD
        % LFM model used for chl prediction was trained on5-200 m depth interval only
        % --> Cant be used for any other

        % % recover Blind coeffs
        % KdFDcoefs = getcoef(fdKd_iP);
        % lumBlindForPred_fd = fd(KdFDcoefs,basisObj);
        %
        % % compute Blind prediction
        % % CENTER THE PREDICTOR DATA
        % xcBlind_fd = center(lumBlindForPred_fd);
        % % COMPUTE PREDICTION
        % yhatpenBlind_fdCoefs_temp = defVals.LFM.LFMcore.Bpen*getcoef(xcBlind_fd);
        %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE CHLA PREDICTION
        % % eval CHLAFITBlind / CHLAHATpBlind
        % newObj_dfCoefs = yhatpenBlind_fdCoefs_temp + defVals.LFM.LFMcore.meanFittedObsFDcoefs;
        % newObj_fd = fd(newObj_dfCoefs, basisObj);
        % platform_processed.CHL_LFM.Full(fitInterval,iP) = eval_fd(newObj_fd,fitInterval);
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREDICTION POST-PROCESSING
        % % remove negative values
        % platform_processed.CHL_LFM.Full(platform_processed.CHL_LFM.Full < 0) = NaN;

    end

    % Compute linear PAR array
    tagProcessed.PAR.lin_RegDrkSatFitAll = exp(tagProcessed.PAR.log_RegDrkSatFitAll);

    %% Functional fit (predefined depth interval)
    tagProcessed.PAR.log_RegDrkSatFitBnd = NaN(size(tagProcessed.PAR.log_RegDrkSat));
    tagProcessed.PAR.Kd_bnd = NaN(size(tagProcessed.PAR.log_RegDrkSat));
    tagProcessed.CHL_LFM.Bound = NaN(size(tagProcessed.PAR.log_RegDrkSat));

    % Select profiles for fit computation
    lumToFit_bnd = tagProcessed.PAR.log_RegDrkSat;
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
        tagProcessed.PAR.log_RegDrkSatFitBnd(fitInterval,i_profiles) = eval_fd(fdFit_PAR,fitInterval);

        % Calculate Kd as the derivative function of the PAR fit
        fdKd = deriv_fd(fdFit_PAR);
        tagProcessed.PAR.Kd_bnd(fitInterval,i_profiles) = -eval_fd(fdKd,fitInterval);

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
        tagProcessed.CHL_LFM.Bound(fitInterval,i_profiles) = eval_fd(newObj_fd,fitInterval);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREDICTION POST-PROCESSING
        % remove negative values
        tagProcessed.CHL_LFM.Bound(tagProcessed.CHL_LFM.Bound < 0) = NaN;
        % fill upper part of profile
        % CHLAupper = repmat(platform_processed.CHL_LFM.Bound(topBound_lfm,:),topBound_lfm - 1,1);
        % CHLA_LFM(1:topBound_lfm - 1,:) = CHLAupper;
    end

    % Compute linear PAR array
    tagProcessed.PAR.lin_RegDrkSatFitBnd = exp(tagProcessed.PAR.log_RegDrkSatFitBnd);

    %% Final things
    % - New subsurface PAR value based on fitted data
    % - Euphotic depth (1% of the surface PAR depth, after Morel and Berthon, 1989)
    iSurfVal_fit = find_ndim(isfinite(tagProcessed.PAR.lin_RegDrkSatFitAll),1,'first');
    finiteIndices = find(iSurfVal_fit>0);
    for iP = finiteIndices
        ProfileInfo_PAR.surfaceValueFit(iP,1) = tagProcessed.PAR.lin_RegDrkSatFitAll(iSurfVal_fit(iP),iP);
        Zeu = -find(tagProcessed.PAR.lin_RegDrkSatFitAll(:,iP) / ProfileInfo_PAR.surfaceValueFit(iP) < 0.01,1,'first');
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
    iZ_quench = find_ndim(tagProcessed.PAR.lin_RegDrkSatFitAll > 15, 1, 'last')';
    iZ_quench(iZ_quench == 0) = NaN;
    ProfileInfo_PAR.quenchDepth = -iZ_quench;

    % (ML-)mean Kd
    ProfileInfo_PAR.meanKd = mean(tagProcessed.PAR.Kd_all,1,'omitnan')';
    for iP = 1 : tagMetadata.nProfs
        iZ_MLD = abs(round(ProfileInfo.MLD(iP)));
        if isfinite(iZ_MLD)
            ProfileInfo_PAR.meanKdML(iP) = mean(tagProcessed.PAR.Kd_all(1:iZ_MLD,iP),1,'omitnan');
        end
    end
end

%% CMD message: done
fprintf('\b\b \x2713\n')

end