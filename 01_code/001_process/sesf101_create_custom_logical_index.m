%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% create custom logicals to filter profiles according to project requirements

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x
% sesf01x
% sesf03x
% sesf04x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
disp('***MATLAB version information***')
disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% last modified: 23.03.31
%
% 20.07.01 update: add idx000_genProfileSelec and idx013_bathyMin to
% custExport
% 21.01.12 update: adjusting variables names (CHLA_nad and chlaData)
% 21.01.28 update: modifying variables names according to new nomenclature
% 22.02.03 updates:
% adding writing of idx variables in genData table
% modifying computation of idx101 & idx102 with lfm-specific top/botBound
% 22.02.05 update:
% removing idx001 (loc) from idx101 & idx102 computation
% 
% -----------------------------------------------------------------------

disp(strcat('creating custom logical indexing for filtering purposes for platform:',...
    platform_metadata.platform_name,'.....'))

%% logical index for profiles to be used in

% LINEAR FUNCTIONAL MODEL TRAINING
idx101_lfmTraining =...
    idx000_genProfileSelec &...
    idx013_bathyMin &...
    idx031_lightNonNan &...
    parData.saturationDepth <= topBound_lfm &...
    parData.darkDepth > botBound_lfm &...
    idx035_lightDay &...
    idx041_chlaNonNan &...
    fluoData.lastNonNan >= botBound_lfm ;
np101_lfmTraining = nnz(idx101_lfmTraining) ;
disp(['dataset:' int2str(np101_lfmTraining)...
    ' profiles to be exported for model training'])

% CHLALUM PREDICTION
idx102_chlalumPredExport =...
    idx000_genProfileSelec &...
    idx013_bathyMin &...
    idx031_lightNonNan &...
    parData.saturationDepth <= topBound_lfm &...
    parData.darkDepth > botBound_lfm &...
    idx035_lightDay ;
np102_chlalumPredExport = nnz(idx102_chlalumPredExport) ;
disp(['dataset:' int2str(np102_chlalumPredExport)...
    ' profiles to be exported for CHLALUM prediction'])

% ANY OTHER REQUIRED EXPORT
idx103_custExport =...
    idx000_genProfileSelec &...
    idx001_latlonNonNan &...
    idx013_bathyMin &...
    idx031_lightNonNan &...
    idx032_lightInterval &...
    idx035_lightDay &...
    idx041_chlaNonNan &...
    idx042_chlaInterval ;
np103_custExport = nnz(idx103_custExport) ;
disp(['dataset:' int2str(np103_custExport)...
    ' profiles to be exported for custom export'])


%% write a few idx variables in genData table

genData.idx000_genProfileSelec = idx000_genProfileSelec ;
genData.idx031_lightNonNan = idx031_lightNonNan ;
genData.idx035_lightDay = idx035_lightDay ;
genData.idx041_chlaNonNan = idx041_chlaNonNan ;
genData.idx100_allParamNonNan = idx100_allParamNonNan ;
genData.idx101_lfmTraining = idx101_lfmTraining ;
genData.idx102_chlalumPredExport = idx102_chlalumPredExport ;

%% clear temp data
clear ans

%% END