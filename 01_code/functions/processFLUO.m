function [Data,ProfileInfo] = processFLUO(Data,ProfileInfo,defaultPars)
% PROCESSPAR
%
% INPUT ARGUMENTS
%
% OUTPUT

%% CMD message: start
fprintf('Processing <strong>fluorescence</strong> data...');

%% Create fluoData info table
var_names = {...
    'Profile', ...              % profile number
    'noData',...                % no data in profile
    'darkDepth', ...            % depth of first value belonging to dark signal
    'darkValue', ...            % dark value (median of 190 m-200 m values)
    'darkValCorr', ...          % dark value for correction (movmedian + linear fit)
    'relativeSD_ML', ...        % relative SD within mixing layer
    'NPQDepth', ...             % depth from which NPQ correction is applied (see sesf045 script)
    'surfaceVal', ...           % surface fluorescence value
    'maxFluoDepth', ...         % depth of fluo max value
    'maxFluoValue', ...         % fluo max value
    'FirstOD_mean',...          % mean fluorescence in the first optical depth
    'ML_mean',...               % mean fluorescence in the mixed layer
    'MLDbio'...                 % MLDbio (Lacour et al. 2017)
    }';
ProfileInfo.FLUO = array2table(NaN(Data.MetaData.nProfs, numel(var_names)),'VariableNames',var_names);

%% Populate initial information
% Profile #
ProfileInfo.FLUO.Profile = (1 : Data.MetaData.nProfs)';

% [DISABLED] Remove values below pre-defined maximum depth
% Data.Processed.FLUO.Reg(defaultPars.depthInterpGrid < defaultPars.CHL.maxCHLA_depth,:) = NaN;

% Set profiles with constant fluorescence values NaN (due to sensor issues)
Data.Processed.FLUO.Reg(:,range(Data.Processed.FLUO.Reg) == 0) = NaN;

% Empty profiles (all NaN or constant values)
ProfileInfo.FLUO.noData = all(isnan(Data.Processed.FLUO.Reg))';

% Set profiles with too few non-zero values NaN (for functional fit analysis, n < nBasis)
n_nonzero = arrayfun(@(a) nnz(~isnan(Data.Processed.FLUO.Reg(:,a))),1 : Data.MetaData.nProfs);
Data.Processed.FLUO.Reg(:,n_nonzero <= defaultPars.LFM.nBasis) = NaN;

% Shallowest and deepest available observation
% Note: Matching non-NaN indices (i_first/last...) to depths in case interpolation depth grid was
% chosen at increments other than 1 m, in which case the indices are not equal to depth.
firstObs = find_ndim(isfinite(Data.Processed.FLUO.Reg),1,'first')';
lastObs = find_ndim(isfinite(Data.Processed.FLUO.Reg),1,'last')';

%% Processing: Dark correction, NPQ correction, spline fitting
% Initiate processed data fields
Data.Processed.FLUO.RegDrk = NaN(size(Data.Processed.DEPTH));
Data.Processed.FLUO.RegDrkNPQ = NaN(size(Data.Processed.DEPTH));
Data.Processed.FLUO.RegDrkNPQFitAll = NaN(size(Data.Processed.DEPTH));
Data.Processed.FLUO.RegDrkNPQFitBnd = NaN(size(Data.Processed.DEPTH));

% Only proceed with processing if there is any usable data
if ~all(ProfileInfo.FLUO.noData)
    %% Dark correction
    % % Dark value: median of 10 deepest fluo values (interpolated to profiles where no dark estimate could be retrieved)
    % darkDepthDelta = 10;
    % fluoDeep = Data.Processed.FLUO.Reg(-defaultPars.CHL.maxCHLA_depth-(darkDepthDelta-1) : -defaultPars.CHL.maxCHLA_depth,:);
    % ProfileInfo.FLUO.darkValue = fillmissing(median(fluoDeep,1,'omitnan'),'nearest')';
    % 
    % if ~all(isnan(ProfileInfo.FLUO.darkValue))
    %     % Dark value smoothing
    %     nDays = 10;   % time window to be used for smoothing (days), start with 10 days and reduce if obs time period is too short
    %     darkValCorr = movmedian(ProfileInfo.FLUO.darkValue,ceil(Data.MetaData.nProfs/max(ProfileInfo.General.DeployDay)*nDays),'Endpoints','fill');
    %     darkValCorr = fillmissing(darkValCorr,'nearest');
    %     while all(isnan(darkValCorr))
    %         % if dark values are all NaN, the median smoothing window was too big
    %         % reduce window by 1 day until dark values are retrieved
    %         nDays = nDays-1;
    %         darkValCorr = movmedian(ProfileInfo.FLUO.darkValue,ceil(Data.MetaData.nProfs/max(ProfileInfo.General.DeployDay)*nDays),'Endpoints','fill');
    %         darkValCorr = fillmissing(darkValCorr,'nearest');
    %     end
    % 
    %     % Approximate drift over time by calculating a linear regression through smoothed dark values
    %     darkValCorrFit = fitlm(datenum(ProfileInfo.General.Date),darkValCorr);
    %     ProfileInfo.FLUO.darkValCorr = darkValCorrFit.Fitted;
    % 
    %     % Remove dark value and set negative values 0
    %     Data.Processed.FLUO.RegDrk = Data.Processed.FLUO.Reg;
    %     fluoOffset = repmat(ProfileInfo.FLUO.darkValCorr',size(Data.Processed.FLUO.RegDrk,1),1);
    %     Data.Processed.FLUO.RegDrk = Data.Processed.FLUO.RegDrk - fluoOffset;
    %     Data.Processed.FLUO.RegDrk(Data.Processed.FLUO.RegDrk<0) = 0;
    % end

    %% Updated dark correction (consistent with PAR, random subsamples)
    % Dark values (median of the 10 deepest observations, or as specified by darkDepthDelta) are determined for a random
    % subsample of deep profiles (>100 m, N = 10 or as defined in defaultPars.PAR.nDarkProfiles) and used to calculate a
    % tag-specific dark value. The randomised subsampling of deep profiles is repeated nRep times to yield a more robust
    % estimate. The tag-specific dark value is then calculated as the median of the median of all individual dark values (per
    % subsample). NB: The random number generator algorithm is reset within the script to ensure that the results do not
    % change when reprocessing data.

    rng(0,'twister')                                % Reset random number generator algorithm (for reproducibility of the random subsampling)
    deepProfs = find(lastObs > 200)';               % Indices of all deep profiles (deeper than 100 m)
    darkDepthDelta = 10;                            % Depth range over which dark depth is calculated (m)
    nRep = 10;                                      % Number of repetitions of random deep profile subsamples
    darkValues = NaN(Data.MetaData.nProfs,nRep);    % Dark value matrix
    for iRep = 1 : nRep
        % Randomly subsample deep profiles to be used for the dark value calculation (N defined in setDefaults)
        nSamples = min(defaultPars.PAR.nDarkProfiles,numel(deepProfs)); % N subsamples cannot be less than the number of deep profiles to sample from
        % iDarkProfiles = datasample(deepProfs,nSamples,'Replace',false); % Indices of subsampled profiles
        iDarkProfiles = datasample(deepProfs,nSamples,'Replace',false); % Indices of subsampled profiles

        darkValues(iDarkProfiles,iRep) = median(Data.Processed.FLUO.Reg(-defaultPars.CHL.maxCHLA_depth-(darkDepthDelta-1) : -defaultPars.CHL.maxCHLA_depth,iDarkProfiles));
    end

    darkValue = median(median(darkValues,1,'omitnan'),2,'omitnan');

    % Remove dark value and set negative values 0, store dark value in ProfileInfo
    Data.Processed.FLUO.RegDrk = Data.Processed.FLUO.Reg - darkValue;
    Data.Processed.FLUO.RegDrk(Data.Processed.FLUO.RegDrk<=0) = NaN;
    Data.Processed.FLUO.RegDrk = fillmissing(Data.Processed.FLUO.RegDrk,'linear',1,'EndValues','none');
    ProfileInfo.FLUO.darkValue(:) = darkValue;

    %% Detect and remove constant part of the FLUO profile at depth
    for iP = 1 : Data.MetaData.nProfs
        for iZ = firstObs(iP) : lastObs(iP) - 1
            % proceed to cst test (dark offset has been previously removed)
            iZ_FLUO_range = range(Data.Processed.FLUO.RegDrk(iZ:lastObs(iP,1),iP));

            if iZ_FLUO_range <= defaultPars.CHL.delta_for_cstCHLA

                Z_const = iZ + 1;  % depth index first value belonging to constant signal
                Data.Processed.FLUO.RegDrk(Z_const:end,iP) = NaN;
                Data.Processed.FLUO.RegDrk(Z_const:end,iP) = NaN;
                lastObs(iP) = Z_const-1;
                ProfileInfo.FLUO.darkDepth(iP) = -Z_const;

                break   % end for loop when FLUO range falls below threshold value where FLUO is considered constant
            end
        end
    end

    %% Compute relative standard deviation
    for iP = 1 : Data.MetaData.nProfs
        if isfinite(ProfileInfo.General.MLD(iP))
            % Calculate FLUO standard deviation between the MLD and the surface or the quenching depth
            botInd = abs(round(ProfileInfo.General.MLD(iP)));

            if isnan(ProfileInfo.PAR.quenchDepth(iP))
                % If the quenching depth is NaN, no quenching was detected (night time profile)
                % calculate the SD/mean begtween the MLD and the surface.
                topInd = 1;

                stdFluo = std(Data.Processed.FLUO.RegDrk(topInd:botInd,iP),0,'omitnan');
                meanFluo = mean(Data.Processed.FLUO.RegDrk(topInd:botInd,iP),'omitnan');
                ProfileInfo.FLUO.relativeSD_ML(iP) = stdFluo / meanFluo;

            elseif abs(round(ProfileInfo.General.MLD(iP))) > abs(round(ProfileInfo.PAR.quenchDepth(iP)))
                % If the MLD is deeper than the quenching depth
                % calculate the SD/mean between the MLD and the quenching depth
                topInd = abs(round(ProfileInfo.PAR.quenchDepth(iP,1)));

                stdFluo = std(Data.Processed.FLUO.RegDrk(topInd:botInd,iP),0,'omitnan');
                meanFluo = mean(Data.Processed.FLUO.RegDrk(topInd:botInd,iP),'omitnan');
                ProfileInfo.FLUO.relativeSD_ML(iP) = stdFluo / meanFluo;

            else
                % If the MLD is shallower than the quenching depth
                % do not calculate the realtive SD.
                ProfileInfo.FLUO.relativeSD_ML(iP) = NaN;
            end

        end
    end

    %% NPQ correction (Xing et al. 2018: X12+ algorithm)
    Data.Processed.FLUO.RegDrkNPQ = Data.Processed.FLUO.RegDrk;

    % Smoothing FLUO
    w = 1; % no smoothing! [\11-point median filter]
    finiteVals = isfinite(Data.Processed.FLUO.RegDrk);
    FLUOsmooth = movmedian(Data.Processed.FLUO.RegDrk,w,'EndPoints','fill');
    FLUOsmooth(finiteVals) = fillmissing(FLUOsmooth(finiteVals),'nearest');

    for iP = 1 : Data.MetaData.nProfs
        % Find the shallower of the MLD*0.9 and quenching depth (>15 mol quanta): "NPQ layer"
        if isfinite(ProfileInfo.PAR.quenchDepth(iP)) || isfinite(ProfileInfo.General.MLD(iP))
            iZ_NPQ = abs(round(max([ProfileInfo.General.MLD(iP)*.9,ProfileInfo.PAR.quenchDepth(iP)])));
        else
            iZ_NPQ = NaN;
        end

        % Only correct if profile was taken during daytime (and if NPQ depth is finite and not 0)
        if ProfileInfo.General.Daytime(iP) && isfinite(iZ_NPQ) && iZ_NPQ > 0
            % FCHLA maximum within the NPQ layer
            [MaxFluo_X12,zMaxFluo_X12] = max(FLUOsmooth(1:iZ_NPQ,iP));
            % correct FLUO_nadDk
            Data.Processed.FLUO.RegDrkNPQ(1:zMaxFluo_X12,iP) = MaxFluo_X12;
            % NPQ layer depth
            ProfileInfo.FLUO.NPQDepth(iP) = -zMaxFluo_X12;
        end
    end

    %% Functional fit (full profile)
    Data.Processed.FLUO.RegDrkNPQFitAll = NaN(size(Data.Processed.FLUO.RegDrkNPQ));

    % Select profiles for fit computation
    fluoToFit_all = Data.Processed.FLUO.RegDrkNPQ;
    i_profiles = find(...
        ProfileInfo.General.Processed &...      % Passed processing checks in loadData
        ~ProfileInfo.FLUO.noData)';           % has data

    % Save all fd coefficients
    fluo_fdAllCoefs = NaN(defaultPars.LFM.nBasis,Data.MetaData.nProfs);

    % Compute fit
    for iP = i_profiles
        % PAR data to compute fit over
        FLUO_fit = fluoToFit_all(:,iP);
        fitInterval = find(isfinite(FLUO_fit));
        FLUO_fit = FLUO_fit(fitInterval);

        if numel(FLUO_fit) >= 30
            % Bspline functional fit
            basisObj = create_bspline_basis([fitInterval(1) fitInterval(end)],defaultPars.LFM.nBasis,defaultPars.LFM.nOrder);
            fdObj = data2fd(fitInterval,FLUO_fit,basisObj);

            % Evaluate fit and store in FLUO structure
            Data.Processed.FLUO.RegDrkNPQFitAll(fitInterval,iP) = eval_fd(fdObj,fitInterval);

            % Store fd coefficients
            fluo_fdAllCoefs(:,iP) = getcoef(fdObj);
        end
    end

    %% Functional fit (predefined depth interval)
    Data.Processed.FLUO.RegDrkNPQFitBnd = NaN(size(Data.Processed.FLUO.RegDrkNPQ));

    % Select profiles for fit computation
    fluoToFit_bnd = Data.Processed.FLUO.RegDrkNPQ;
    i_profiles = find(...
        (firstObs <= abs(defaultPars.upperDepthBound)... % Finite data within bounds
        & lastObs >= abs(defaultPars.lowerDepthBound)) &...
        ProfileInfo.General.Processed);            % Passed processing checks in loadData

    % proceed only if any profiles were found to have finite PAR data within the specified bounds
    if ~isempty(i_profiles)
        fitInterval = (abs(defaultPars.upperDepthBound) : abs(defaultPars.lowerDepthBound))';
        FLUO_fit = fluoToFit_bnd(fitInterval,i_profiles);

        % Bspline functional fit
        basisObj = create_bspline_basis([fitInterval(1) fitInterval(end)],defaultPars.LFM.nBasis,defaultPars.LFM.nOrder);
        fdObj = data2fd(fitInterval,FLUO_fit,basisObj);

        % Evaluate fit and store in FLUO structure
        Data.Processed.FLUO.RegDrkNPQFitBnd(fitInterval,i_profiles) = eval_fd(fdObj,fitInterval);
    end

    %% Final things
    % Surface fluorescence value
    for iP = 1 : Data.MetaData.nProfs
        firstObs = find_ndim(isfinite(Data.Processed.FLUO.RegDrkNPQ(:,iP)),1,'first');
        if firstObs > 0
            ProfileInfo.FLUO.surfaceVal(iP) = Data.Processed.FLUO.RegDrkNPQ(firstObs,iP)';
        end
    end

    % Mean fluorescence over the first optical depth (for comparison w/ satellite obs)
    iZ_1stOD = abs(round(ProfileInfo.PAR.FirstOD));
    ProfileInfo.FLUO.FirstOD_mean(find(isfinite(iZ_1stOD))) =...
        arrayfun(@(a) mean(Data.Processed.FLUO.RegDrkNPQ(1:iZ_1stOD(a),a),'omitnan'),...
        find(isfinite(ProfileInfo.PAR.FirstOD)))';

    % Mean fluorescence over the mixed layer
    iZ_MLD = abs(round(ProfileInfo.General.MLD));
    ProfileInfo.FLUO.ML_mean(find(isfinite(ProfileInfo.General.MLD))) =...
        arrayfun(@(a) mean(Data.Processed.FLUO.RegDrkNPQ(1:iZ_MLD(a),a),'omitnan'),...
        find(isfinite(ProfileInfo.General.MLD)))';

    % Fluorscence maximum and depth
    [maxFluoValue,Z_maxFluo] = max(Data.Processed.FLUO.RegDrkNPQFitAll);
    ProfileInfo.FLUO.maxFluoValue = maxFluoValue';
    ProfileInfo.FLUO.maxFluoDepth(isfinite(maxFluoValue)) = -Z_maxFluo(isfinite(maxFluoValue))';

    % MLDbio
    i_profiles = find(...
        ProfileInfo.General.Processed &...              % Passed processing checks in loadData
        sum(~isnan(Data.Processed.FLUO.RegDrkNPQFitAll))' >= 3)';  % At least three values

    for iP = i_profiles
        fdCoeffs = fluo_fdAllCoefs(:,iP);
        fitInterval = find(~isnan(Data.Processed.FLUO.RegDrkNPQFitAll(:,iP)));
        basisObj = create_bspline_basis([fitInterval(1) fitInterval(end)],defaultPars.LFM.nBasis,defaultPars.LFM.nOrder);
        fdObj = fd(fdCoeffs,basisObj);
        fdObjDeriv = deriv_fd(fdObj);
        FluoDerivative = eval_fd(fdObjDeriv,fitInterval);
        [~,iMinFDeriv] = min(FluoDerivative);
        ProfileInfo.FLUO.MLDbio(iP) = -iMinFDeriv;
    end
end

%% CMD message: done
fprintf('\b\b \x2713\n\n')

end