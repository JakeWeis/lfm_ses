%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% predict Chl-a from Kd profiles

%%%%%%%%%% REFERENCES %%%%%%%%%%

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% J. Weis (jw)
% last modified: 22.10.28
% 
% ----------------------------------------------------------------------

disp(strcat('running LFM model fine scale PREDICTION'))


%% run model with nbaz ranging in nbazRange

mybL_temp = LFMcore.mybL ;

% recover Blind coeffs
lumBlindForPred_fd = fd(KdFDcoefs(:,idx036_lightFDAfitBnd),...
    mybL_temp) ;
    
    
% compute Blind prediction
% CENTER THE PREDICTOR DATA
xcBlind_fd = center(lumBlindForPred_fd) ;
% PREDICTOR MEAN
xmBlind_fd = mean(lumBlindForPred_fd) ;
% COMPUTE PREDICTION
yhatpenBlind_fdCoefs_temp = LFMcore.Bpen*getcoef(xcBlind_fd) ;   

% eval CHLAFITBlind / CHLAHATpBlind
% INITIALIZE DATA
CHLA_LFM = nan(max_depth,platform_metadata.np) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE CHLA PREDICTION
newObj_dfCoefs = yhatpenBlind_fdCoefs_temp + LFMcore.meanFittedObsFDcoefs ;
newObj_fd = fd(newObj_dfCoefs, mybL_temp) ;
CHLA_LFM(pp,idx036_lightFDAfitBnd) = eval_fd(newObj_fd,pp) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREDICTION POST-PROCESSING
% remove negative values
indNeg_temp = CHLA_LFM < 0 ;
CHLA_LFM(indNeg_temp) = NaN ;
% fill upper part of profile
CHLAupper_temp = repmat(CHLA_LFM(topBound_lfm,:),topBound_lfm - 1,1) ;
CHLA_LFM(1:topBound_lfm - 1,:) = CHLAupper_temp ;


%% clear temp variables
clear *_temp


%% END

