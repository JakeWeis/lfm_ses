function [tagProcessed,ProfileInfo_FLUO] = processFLUO(tagMetadata,tagProcessed,ProfileInfo,ProfileInfo_PAR,defaultPars)
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
    'NPQDepth', ...         % depth from which NPQ correction is applied (see sesf045 script)
    'surfaceVal', ...           % surface fluorescence value
    'maxFluoDepth', ...         % depth of fluo max value
    'maxFluoValue', ...         % fluo max value
    'FirstOD_mean',...          % mean fluorescence in the first optical depth
    'ML_mean',...               % mean fluorescence in the mixed layer
    'MLDbio'...                 % MLDbio (Lacour et al. 2017)
    }';
ProfileInfo_FLUO = array2table(NaN(tagMetadata.nProfs, numel(var_names)),'VariableNames',var_names);

%% Populate initial information
% Profile #
ProfileInfo_FLUO.Profile = (1 : tagMetadata.nProfs)';

% Remove values below pre-defined maximum depth
tagProcessed.FLUO.Reg(defaultPars.depthInterpGrid < defaultPars.CHL.maxCHLA_depth,:) = NaN;

% Set profiles with constant fluorescence values NaN (due to sensor issues)
tagProcessed.FLUO.Reg(:,range(tagProcessed.FLUO.Reg) == 0) = NaN;

% Empty profiles (all NaN or constant values)
ProfileInfo_FLUO.noData = all(isnan(tagProcessed.FLUO.Reg))';

% Set profiles with too few non-zero values NaN (for functional fit analysis, n < nBasis)
n_nonzero = arrayfun(@(a) nnz(~isnan(tagProcessed.FLUO.Reg(:,a))),1 : tagMetadata.nProfs);
tagProcessed.FLUO.Reg(:,n_nonzero <= defaultPars.LFM.nBasis) = NaN;

% Shallowest and deepest available observation
% Note: Matching non-NaN indices (i_first/last...) to depths in case interpolation depth grid was
% chosen at increments other than 1 m, in which case the indices are not equal to depth.
firstObs = find_ndim(isfinite(tagProcessed.FLUO.Reg),1,'first')';
lastObs = find_ndim(isfinite(tagProcessed.FLUO.Reg),1,'last')';

% Only proceed with processing if there is any usable data
if ~all(ProfileInfo_FLUO.noData)
    %% Dark correction
    % Dark value: median of 10 deepest fluo values (interpolated to profiles where no dark estimate could be retrieved)
    darkDepthDelta = 10;
    fluoDeep = tagProcessed.FLUO.Reg(-defaultPars.CHL.maxCHLA_depth-(darkDepthDelta-1) : -defaultPars.CHL.maxCHLA_depth,:);
    ProfileInfo_FLUO.darkValue = fillmissing(median(fluoDeep,1,'omitnan'),'nearest')';

    % Dark value smoothing
    nDays = 10;   % time window to be used for smoothing (days), start with 10 days and reduce if obs time period is too short
    darkValCorr = movmedian(ProfileInfo_FLUO.darkValue,ceil(tagMetadata.nProfs/max(ProfileInfo.DeployDay)*nDays),'Endpoints','fill');
    darkValCorr = fillmissing(darkValCorr,'nearest');
    while all(isnan(darkValCorr))
        % if dark values are all NaN, the median smoothing window was too big
        % reduce window by 1 day until dark values are retrieved
        nDays = nDays-1;
        darkValCorr = movmedian(ProfileInfo_FLUO.darkValue,ceil(tagMetadata.nProfs/max(ProfileInfo.DeployDay)*nDays),'Endpoints','fill');
        darkValCorr = fillmissing(darkValCorr,'nearest');
    end

    % Approximate drift over time by calculating a linear regression through smoothed dark values
    darkValCorrFit = fitlm(datenum(ProfileInfo.Date),darkValCorr);
    ProfileInfo_FLUO.darkValCorr = darkValCorrFit.Fitted;

    % Remove dark value and set negative values 0
    tagProcessed.FLUO.RegDrk = tagProcessed.FLUO.Reg;
    fluoOffset = repmat(ProfileInfo_FLUO.darkValCorr',size(tagProcessed.FLUO.RegDrk,1),1);
    tagProcessed.FLUO.RegDrk = tagProcessed.FLUO.RegDrk - fluoOffset;
    tagProcessed.FLUO.RegDrk(tagProcessed.FLUO.RegDrk<0) = 0;


    %% Detect and remove constant part of the FLUO profile at depth
    for iP = 1 : tagMetadata.nProfs
        for iZ = firstObs(iP) : lastObs(iP) - 1
            % proceed to cst test (dark offset has been previously removed)
            iZ_FLUO_range = range(tagProcessed.FLUO.RegDrk(iZ:lastObs(iP,1),iP));

            if iZ_FLUO_range <= defaultPars.CHL.delta_for_cstCHLA

                Z_const = iZ + 1;  % depth index first value belonging to constant signal
                tagProcessed.FLUO.RegDrk(Z_const:end,iP) = NaN;
                tagProcessed.FLUO.RegDrk(Z_const:end,iP) = NaN;
                lastObs(iP) = Z_const-1;
                ProfileInfo_FLUO.darkDepth(iP) = -Z_const;

                break   % end for loop when FLUO range falls below threshold value where FLUO is considered constant
            end
        end
    end

    %% Compute relative standard deviation
    for iP = 1 : tagMetadata.nProfs
        if isfinite(ProfileInfo.MLD(iP))
            % Calculate FLUO standard deviation between the MLD and the surface or the quenching depth
            botInd = abs(round(ProfileInfo.MLD(iP)));

            if isnan(ProfileInfo_PAR.quenchDepth(iP))
                % If the quenching depth is NaN, no quenching was detected (night time profile)
                % calculate the SD/mean begtween the MLD and the surface.
                topInd = 1;

                stdFluo = std(tagProcessed.FLUO.RegDrk(topInd:botInd,iP),0,'omitnan');
                meanFluo = mean(tagProcessed.FLUO.RegDrk(topInd:botInd,iP),'omitnan');
                ProfileInfo_FLUO.relativeSD_ML(iP) = stdFluo / meanFluo;

            elseif abs(round(ProfileInfo.MLD(iP))) > abs(round(ProfileInfo_PAR.quenchDepth(iP)))
                % If the MLD is deeper than the quenching depth
                % calculate the SD/mean between the MLD and the quenching depth
                topInd = abs(round(ProfileInfo_PAR.quenchDepth(iP,1)));

                stdFluo = std(tagProcessed.FLUO.RegDrk(topInd:botInd,iP),0,'omitnan');
                meanFluo = mean(tagProcessed.FLUO.RegDrk(topInd:botInd,iP),'omitnan');
                ProfileInfo_FLUO.relativeSD_ML(iP) = stdFluo / meanFluo;

            else
                % If the MLD is shallower than the quenching depth
                % do not calculate the realtive SD.
                ProfileInfo_FLUO.relativeSD_ML(iP) = NaN;
            end

        end
    end

    %% NPQ correction (Xing et al. 2018: X12+ algorithm)
    tagProcessed.FLUO.RegDrkNPQ = tagProcessed.FLUO.RegDrk;

    % Smoothing FLUO
    w = 1; % no smoothing! [\11-point median filter]
    finiteVals = isfinite(tagProcessed.FLUO.RegDrk);
    FLUOsmooth = movmedian(tagProcessed.FLUO.RegDrk,w,'EndPoints','fill');
    FLUOsmooth(finiteVals) = fillmissing(FLUOsmooth(finiteVals),'nearest');

    for iP = 1 : tagMetadata.nProfs
        % Find the shallower of the MLD and quenching depth (>15 mol quanta): "NPQ layer"
        if isfinite(ProfileInfo_PAR.quenchDepth(iP)) && isfinite(ProfileInfo.MLD(iP))
            iZ_NPQ = abs(round(max([ProfileInfo.MLD(iP),ProfileInfo_PAR.quenchDepth(iP)])));
        else
            iZ_NPQ = NaN;
        end

        % Only correct if profile was taken during daytime (and if NPQ depth is finite and not 0)
        if ProfileInfo.Daytime(iP) && isfinite(iZ_NPQ) && iZ_NPQ > 0
            % FCHLA maximum within the NPQ layer
            [MaxFluo_X12,zMaxFluo_X12] = max(FLUOsmooth(1:iZ_NPQ,iP));
            % correct FLUO_nadDk
            tagProcessed.FLUO.RegDrkNPQ(1:zMaxFluo_X12,iP) = MaxFluo_X12;
            % NPQ layer depth
            ProfileInfo_FLUO.NPQDepth(iP) = -zMaxFluo_X12;
        end
    end

    %% Functional fit (full profile)
    tagProcessed.FLUO.RegDrkNPQFitAll = NaN(size(tagProcessed.FLUO.RegDrkNPQ));

    % Select profiles for fit computation
    fluoToFit_all = tagProcessed.FLUO.RegDrkNPQ;
    i_profiles = find(...
        ProfileInfo.Processed &...      % Passed processing checks in loadData
        ~ProfileInfo_FLUO.noData)';           % has data

    % Save all fd coefficients
    fluo_fdAllCoefs = NaN(defaultPars.LFM.nBasis,tagMetadata.nProfs);

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
            tagProcessed.FLUO.RegDrkNPQFitAll(fitInterval,iP) = eval_fd(fdObj,fitInterval);

            % Store fd coefficients
            fluo_fdAllCoefs(:,iP) = getcoef(fdObj);
        end
    end

    %% Functional fit (predefined depth interval)
    tagProcessed.FLUO.RegDrkNPQFitBnd = NaN(size(tagProcessed.FLUO.RegDrkNPQ));

    % Select profiles for fit computation
    fluoToFit_bnd = tagProcessed.FLUO.RegDrkNPQ;
    i_profiles = find(...
        (firstObs <= abs(defaultPars.upperDepthBound)... % Finite data within bounds
        & lastObs >= abs(defaultPars.lowerDepthBound)) &...
        ProfileInfo.Processed);            % Passed processing checks in loadData

    % proceed only if any profiles were found to have finite PAR data within the specified bounds
    if ~isempty(i_profiles)
        fitInterval = (abs(defaultPars.upperDepthBound) : abs(defaultPars.lowerDepthBound))';
        FLUO_fit = fluoToFit_bnd(fitInterval,i_profiles);

        % Bspline functional fit
        basisObj = create_bspline_basis([fitInterval(1) fitInterval(end)],defaultPars.LFM.nBasis,defaultPars.LFM.nOrder);
        fdObj = data2fd(fitInterval,FLUO_fit,basisObj);

        % Evaluate fit and store in FLUO structure
        tagProcessed.FLUO.RegDrkNPQFitBnd(fitInterval,i_profiles) = eval_fd(fdObj,fitInterval);
    end

    %% Final things
    % Surface fluorescence value
    for iP = 1 : tagMetadata.nProfs
        firstObs = find_ndim(isfinite(tagProcessed.FLUO.RegDrkNPQ(:,iP)),1,'first');
        if firstObs > 0
            ProfileInfo_FLUO.surfaceVal(iP) = tagProcessed.FLUO.RegDrkNPQ(firstObs,iP)';
        end
    end

    % Mean fluorescence over the first optical depth (for comparison w/ satellite obs)
    iZ_1stOD = abs(round(ProfileInfo_PAR.FirstOD));
    ProfileInfo_FLUO.FirstOD_mean(find(isfinite(iZ_1stOD))) =...
        arrayfun(@(a) mean(tagProcessed.FLUO.RegDrkNPQ(1:iZ_1stOD(a),a),'omitnan'),...
        find(isfinite(ProfileInfo_PAR.FirstOD)))';

    % Mean fluorescence over the mixed layer
    iZ_MLD = abs(round(ProfileInfo.MLD));
    ProfileInfo_FLUO.ML_mean(find(isfinite(ProfileInfo.MLD))) =...
        arrayfun(@(a) mean(tagProcessed.FLUO.RegDrkNPQ(1:iZ_MLD(a),a),'omitnan'),...
        find(isfinite(ProfileInfo.MLD)))';

    % Fluorscence maximum and depth
    [maxFluoValue,Z_maxFluo] = max(tagProcessed.FLUO.RegDrkNPQFitAll);
    ProfileInfo_FLUO.maxFluoValue = maxFluoValue';
    ProfileInfo_FLUO.maxFluoDepth(isfinite(maxFluoValue)) = -Z_maxFluo(isfinite(maxFluoValue))';

    % MLDbio
    i_profiles = find(...
        ProfileInfo.Processed &...              % Passed processing checks in loadData
        sum(~isnan(tagProcessed.FLUO.RegDrkNPQFitAll))' >= 3)';  % At least three values

    for iP = i_profiles
        fdCoeffs = fluo_fdAllCoefs(:,iP);
        fitInterval = find(~isnan(tagProcessed.FLUO.RegDrkNPQFitAll(:,iP)));
        basisObj = create_bspline_basis([fitInterval(1) fitInterval(end)],defaultPars.LFM.nBasis,defaultPars.LFM.nOrder);
        fdObj = fd(fdCoeffs,basisObj);
        fdObjDeriv = deriv_fd(fdObj);
        FluoDerivative = eval_fd(fdObjDeriv,fitInterval);
        [~,iMinFDeriv] = min(FluoDerivative);
        ProfileInfo_FLUO.MLDbio(iP) = -iMinFDeriv;
    end
end

%% CMD message: done
fprintf('\b\b \x2713\n\n')

end