%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% compute CHLA data raw metrics:
% - definition interval

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x
% sesf01x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
disp('***MATLAB version information***')
disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% last modified: 22.03.29
%
% -----------------------------------------------------------------------

disp(strcat('FLUO data pre-processing step1: raw METRICS / platform:',...
    platform_metadata.platform_name,'.....'))


%% compute definition interval (first/last non NaN values of profile)

% first/last nonNaN value of each profile (based on CHLA_nadReg array)
firstNonNan_temp = find_ndim(~isnan(FLUO_nadReg),1,'first').' ;
firstNonNan_temp(firstNonNan_temp == 0) = NaN ;
fluoData.firstNonNan = firstNonNan_temp ;
lastNonNan_temp = find_ndim(~isnan(FLUO_nadReg),1,'last').' ;
lastNonNan_temp(lastNonNan_temp == 0) = NaN ;
fluoData.lastNonNan = lastNonNan_temp ;


%% clear temp data
clear *_temp


%% END
