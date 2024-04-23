function [platform_processed,fluoData] = processFLUO(platform_metadata,platform_processed,genData,defVals,parData)
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
fluoData = array2table(NaN(platform_metadata.nProfs, numel(var_names)),'VariableNames',var_names);

%% Populate initial information
% Profile #
fluoData.Profile = (1 : platform_metadata.nProfs)';

% Remove values below pre-defined maximum depth
platform_processed.FLUO.Reg(defVals.depthInterpGrid < defVals.CHL.maxCHLA_depth,:) = NaN;

% Empty profiles (all NaN or constant values)
fluoData.noData = all(isnan(platform_processed.FLUO.Reg))' | range(platform_processed.FLUO.Reg,1)' == 0;

% Set profiles with too few non-zero values NaN (for functional fit analysis, n < nBasis)
n_nonzero = arrayfun(@(a) nnz(~isnan(platform_processed.FLUO.Reg(:,a))),1 : platform_metadata.nProfs);
platform_processed.FLUO.Reg(:,n_nonzero <= defVals.LFM.nBasis) = NaN;

% Shallowest and deepest available observation
% Note: Matching non-NaN indices (i_first/last...) to depths in case interpolation depth grid was 
% chosen at increments other than 1 m, in which case the indices are not equal to depth.
firstObs = find_ndim(isfinite(platform_processed.FLUO.Reg),1,'first')';
lastObs = find_ndim(isfinite(platform_processed.FLUO.Reg),1,'last')';


%% Dark correction
% Dark value: median of 10 deepest fluo values (interpolated to profiles where no dark estimate could be retrieved)
darkDepthDelta = 5;
fluoDeep = platform_processed.FLUO.Reg(-defVals.CHL.maxCHLA_depth-(darkDepthDelta-1) : -defVals.CHL.maxCHLA_depth,:);
fluoData.darkValue = fillmissing(median(fluoDeep,1,'omitnan'),'nearest')';

% Dark value smoothing
nDays = 10;     % time window to be used for smoothing (days)
darkValCorr = movmedian(fluoData.darkValue,ceil(platform_metadata.nProfs/max(genData.DeployDay)*nDays),'Endpoints','fill');
darkValCorr = fillmissing(darkValCorr,'nearest');

% Approximate drift over time by calculating a linear regression through smoothed dark values
darkValCorrFit = fitlm(datenum(genData.Date),darkValCorr);
fluoData.darkValCorr = darkValCorrFit.Fitted;

% Remove dark value and set negative values 0
platform_processed.FLUO.RegDrk = platform_processed.FLUO.Reg;
fluoOffset = repmat(fluoData.darkValCorr',size(platform_processed.FLUO.RegDrk,1),1);
platform_processed.FLUO.RegDrk = platform_processed.FLUO.RegDrk - fluoOffset;
platform_processed.FLUO.RegDrk(platform_processed.FLUO.RegDrk<0) = 0;


%% Detect and remove constant part of the FLUO profile at depth
for iP = 1 : platform_metadata.nProfs
    for iZ = firstObs(iP) : lastObs(iP) - 1
        % proceed to cst test (dark offset has been previously removed)
        iZ_FLUO_range = range(platform_processed.FLUO.RegDrk(iZ:lastObs(iP,1),iP));
        
        if iZ_FLUO_range <= defVals.CHL.delta_for_cstCHLA

            Z_const = iZ + 1;  % depth index first value belonging to constant signal
            platform_processed.FLUO.RegDrk(Z_const:end,iP) = NaN;
            platform_processed.FLUO.RegDrk(Z_const:end,iP) = NaN;
            lastObs(iP) = Z_const-1;
            fluoData.darkDepth(iP) = -Z_const;

            break   % end for loop when FLUO range falls below threshold value where FLUO is considered constant
        end
    end
end

%% Compute relative standard deviation
for iP = 1 : platform_metadata.nProfs
    if isfinite(genData.MLD(iP))
        % Calculate FLUO standard deviation between the MLD and the surface or the quenching depth
        botInd = abs(round(genData.MLD(iP)));
        
        if isnan(parData.quenchDepth(iP))
            % If the quenching depth is NaN, no quenching was detected (night time profile)
            % calculate the SD/mean begtween the MLD and the surface.
            topInd = 1;
            
            stdFluo = std(platform_processed.FLUO.RegDrk(topInd:botInd,iP),0,'omitnan');
            meanFluo = mean(platform_processed.FLUO.RegDrk(topInd:botInd,iP),'omitnan');
            fluoData.relativeSD_ML(iP) = stdFluo / meanFluo;

        elseif abs(round(genData.MLD(iP))) > abs(round(parData.quenchDepth(iP)))
            % If the MLD is deeper than the quenching depth
            % calculate the SD/mean between the MLD and the quenching depth
            topInd = abs(round(parData.quenchDepth(iP,1)));
            
            stdFluo = std(platform_processed.FLUO.RegDrk(topInd:botInd,iP),0,'omitnan');
            meanFluo = mean(platform_processed.FLUO.RegDrk(topInd:botInd,iP),'omitnan');
            fluoData.relativeSD_ML(iP) = stdFluo / meanFluo;

        else
            % If the MLD is shallower than the quenching depth
            % do not calculate the realtive SD.
            fluoData.relativeSD_ML(iP) = NaN;
        end
        
    end
end

%% NPQ correction (Xing et al. 2018: X12+ algorithm)
platform_processed.FLUO.RegDrkNPQ = platform_processed.FLUO.RegDrk;

% Smoothing FLUO
w = 11; % 11-point median filter
finiteVals = isfinite(platform_processed.FLUO.RegDrk);
FLUOsmooth = movmedian(platform_processed.FLUO.RegDrk,w,'EndPoints','fill');
FLUOsmooth(finiteVals) = fillmissing(FLUOsmooth(finiteVals),'nearest');

% Find the shallower of the MLD and quenching depth (>15 mol quanta): "NPQ layer"
iZ_NPQ = abs(round(min(genData.MLD(:,1),parData.quenchDepth(:,1))));

for iP = 1 : platform_metadata.nProfs
    % Only correct if profile was taken during daytime
    if genData.Daytime(iP)
        % FCHLA maximum within the NPQ layer
        [MaxFluo_X12,zMaxFluo_X12] = max(FLUOsmooth(1:iZ_NPQ(iP),iP));
        % correct FLUO_nadDk
        platform_processed.FLUO.RegDrkNPQ(1:zMaxFluo_X12,iP) = MaxFluo_X12;
        % NPQ layer depth
        fluoData.NPQDepth(iP) = -zMaxFluo_X12;
    end
end

%% Functional fit (full profile)
platform_processed.FLUO.RegDrkNPQFitAll = NaN(size(platform_processed.FLUO.RegDrkNPQ));

% Select profiles for fit computation
fluoToFit_all = platform_processed.FLUO.RegDrkNPQ;
i_profiles = find(...
    genData.Processed &...      % Passed processing checks in loadData
    ~fluoData.noData)';           % has data

% Save all fd coefficients
fluo_fdAllCoefs = NaN(defVals.LFM.nBasis,platform_metadata.nProfs);

% Compute fit
for iP = i_profiles
    % PAR data to compute fit over
    FLUO_fit = fluoToFit_all(:,iP);
    fitInterval = find(isfinite(FLUO_fit));
    FLUO_fit = FLUO_fit(fitInterval);

    % Bspline functional fit
    basisObj = create_bspline_basis([fitInterval(1) fitInterval(end)],defVals.LFM.nBasis,defVals.LFM.nOrder);
    fdObj = data2fd(fitInterval,FLUO_fit,basisObj);

    % Evaluate fit and store in FLUO structure
    platform_processed.FLUO.RegDrkNPQFitAll(fitInterval,iP) = eval_fd(fdObj,fitInterval);

    % Store fd coefficients
    fluo_fdAllCoefs(:,iP) = getcoef(fdObj);
end

%% Functional fit (predefined depth interval)
platform_processed.FLUO.RegDrkNPQFitBnd = NaN(size(platform_processed.FLUO.RegDrkNPQ));

% Select profiles for fit computation
fluoToFit_bnd = platform_processed.PAR.log_RegDrkSat;
i_profiles = find(...
    (firstObs <= abs(defVals.upperDepthBound)... % Finite data within bounds
    & lastObs >= abs(defVals.lowerDepthBound)) &...
    genData.Processed);            % Passed processing checks in loadData

% proceed only if any profiles were found to have finite PAR data within the specified bounds
if ~isempty(i_profiles)
    fitInterval = (abs(defVals.upperDepthBound) : abs(defVals.lowerDepthBound))';
    FLUO_fit = fluoToFit_bnd(fitInterval,i_profiles);

    % Bspline functional fit
    basisObj = create_bspline_basis([fitInterval(1) fitInterval(end)],defVals.LFM.nBasis,defVals.LFM.nOrder);
    fdObj = data2fd(fitInterval,FLUO_fit,basisObj);

    % Evaluate fit and store in FLUO structure
    platform_processed.FLUO.RegDrkNPQFitBnd(fitInterval,i_profiles) = eval_fd(fdObj,fitInterval);
end

%% Final things
% Surface fluorescence value
for iP = 1 : platform_metadata.nProfs
    firstObs = find_ndim(isfinite(platform_processed.FLUO.RegDrkNPQ(:,iP)),1,'first');
    if firstObs > 0
        fluoData.surfaceVal(iP) = platform_processed.FLUO.RegDrkNPQ(firstObs,iP)';
    end
end

% Mean fluorescence over the first optical depth (for comparison w/ satellite obs)
iZ_1stOD = abs(round(parData.FirstOD));
fluoData.FirstOD_mean(find(isfinite(iZ_1stOD))) =...
    arrayfun(@(a) mean(platform_processed.FLUO.RegDrkNPQ(1:iZ_1stOD(a),a),'omitnan'),...
    find(isfinite(parData.FirstOD)))';

% Mean fluorescence over the mixed layer
iZ_MLD = abs(round(genData.MLD));
fluoData.ML_mean(find(isfinite(genData.MLD))) =...
    arrayfun(@(a) mean(platform_processed.FLUO.RegDrkNPQ(1:iZ_MLD(a),a),'omitnan'),...
    find(isfinite(genData.MLD)))';

% Fluorscence maximum and depth
[maxFluoValue,Z_maxFluo] = max(platform_processed.FLUO.RegDrkNPQFitAll);
fluoData.maxFluoValue = maxFluoValue';
fluoData.maxFluoDepth(isfinite(maxFluoValue)) = -Z_maxFluo(isfinite(maxFluoValue))';

% MLDbio
i_profiles = find(...
    genData.Processed &...              % Passed processing checks in loadData
    sum(~isnan(platform_processed.FLUO.RegDrkNPQFitAll))' >= 3)';  % At least three values

for iP = i_profiles
    fdCoeffs = fluo_fdAllCoefs(:,iP);
    fitInterval = find(~isnan(platform_processed.FLUO.RegDrkNPQFitAll(:,iP)));
    basisObj = create_bspline_basis([fitInterval(1) fitInterval(end)],defVals.LFM.nBasis,defVals.LFM.nOrder);
    fdObj = fd(fdCoeffs,basisObj);
    fdObjDeriv = deriv_fd(fdObj);
    FluoDerivative = eval_fd(fdObjDeriv,fitInterval);
    [~,iMinFDeriv] = min(FluoDerivative);
    fluoData.MLDbio(iP) = -iMinFDeriv;
end

%% CMD message: done
fprintf('\b\b \x2713\n')

end