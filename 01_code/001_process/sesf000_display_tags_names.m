%% script metadata
%
% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% display the names of the tags to use
%
%%%%%%%%%% REFERENCES %%%%%%%%%%
%
%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% none
%
%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% J. Weis (jw)
% last modified: 22/03/2024
% 
% -----------------------------------------------------------------------

%% platform type / data folder / deployment metadata

plfrm_temp = 's' ;
platform_type = 'sealtag' ;
DEPS_METADATA_temp  = readtable([root_data_seal 'DEP_METADATA.csv']) ;         
dataset_type = 'lr' ;      
fold_info = dir(strcat(root_data,'*.nc')) ;

%% Base station

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE HERE THE BASE STATION FOR DATASET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deploymentRef_temp = 'ft22' ; % just for metadata, no impact on processing
% define base station: DDU/PAF/P.Del.
base_station = char(DEPS_METADATA_temp.base_station(strcmp(DEPS_METADATA_temp.deployment_code,...
    deploymentRef_temp) == true)) ;
if isempty(base_station)
    base_station = NaN ;
end

%% display deployment metadata

% display deployment metadata for quick check
deploymentMetadata_temp = table() ;
% deploymentMetadata_temp.deployment_ref = deploymentRef_temp ;
deploymentMetadata_temp.platform_type = platform_type ;
deploymentMetadata_temp.base_station = base_station ;
% deploymentMetadata_temp.year = readtable('DEP_YEAR.csv')
deploymentMetadata_temp.n_platforms = numel(fold_info) ;

deploymentMetadata_temp

displayDeploymentFiles_temp = 'y' ;
deploymentFiles_temp = struct2table(fold_info) ;
if numel(fold_info) <= 1
    deploymentFiles_temp = {[1:numel(fold_info)]' deploymentFiles_temp.name} ;
    deploymentFiles_temp = cell2table(deploymentFiles_temp,...
        'VariableNames',{'num_in_fold','file_name'}) ;
else
deploymentFiles_temp = table([1:numel(fold_info)]',...
    deploymentFiles_temp.name,...
    'VariableNames',{'num_in_fold','file_name'}) ;
end

disp(deploymentFiles_temp)