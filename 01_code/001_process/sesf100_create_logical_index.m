%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% create logicals to filter profiles according to following parameters:
% where all parameters (LAT/LON,FLUO,LIGHT) are non NaN
% where integration interval boundaries are min_topIntegBound<.<min_botIntegBound or larger

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x
% sesf01x
% sesf03x
% sesf04x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% last modified: 22.04.06
% 
% 21.01.12 update: adjusting variables names (CHLA_nad and chlaData)
% 21.01.28 update: modifying variables names according to new nomenclature
% 21.04.07 update: removing chlaData.firstNonNan <= min_topIntegBound
% criteria in idx042 because NPQ correction ensures profile definition from
% depth = 1m

% 
% -----------------------------------------------------------------------

disp(strcat('creating default logical indexing for filtering purposes for platform:',...
    platform_metadata.platform_name,'.....'))


%% create logical indexing where all parameters (LAT/LON,CHLA,LIGHT) are non NaN

% profiles where all parameters (LAT/LON,CHLA,LIGHT) are non NaN
idx100_allParamNonNan =...
    idx001_latlonNonNan &...        % where LAT/LON is available
    idx031_lightNonNan &...            % where LIGHT is non NaN
    idx041_chlaNonNan ;          % where CHLA is non NaN

    
    
np100_allParamNonNan= nnz(idx100_allParamNonNan) ;
disp(strcat('dataset:',...
    int2str(np100_allParamNonNan),...
    ' non NaN (LAT/LON,FLUO,LIGHT) profiles'))




%% create logical indexing for integration interval
% where boundaries are min_topIntegBound<.<min_botIntegBound or larger

% for LIGHT
idx032_lightInterval =...
    parData.saturationDepth <= min_topIntegBound &...
    parData.darkDepth > min_botIntegBound ;
 
np032_lightInterval = nnz(idx031_lightNonNan & idx032_lightInterval) ;
disp(['dataset:' int2str(np032_lightInterval)...
    ' non NaN LIGHT profiles with integration interval larger than ['...
    int2str(min_topIntegBound) ';' int2str(min_botIntegBound) ']'])



% for FLUO
idx042_chlaInterval =...
    fluoData.lastNonNan >= min_botIntegBound ;
% (dark correction does not modify lastnonnan)
    
np042_chlaInterval = nnz(idx041_chlaNonNan & idx042_chlaInterval) ;
disp(['dataset:' int2str(np042_chlaInterval)...
    ' non NaN FLUO profiles with integration interval larger than ['...
    int2str(min_topIntegBound) ';' int2str(min_botIntegBound) ']'])




%% clear temp data
clear ans

%% END