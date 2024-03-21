%% script metadata
%
% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% load platform data
%
%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define platform.m
% sesf001_choose_your_menu
%
%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])
%
%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% J. Weis (jw)
% last modified: 22.02.02
% 
% 20.09.30 updates:
% removing definition of working folders (now in sesf000_define_platform)
% adding quick display of dataset information
% 20.11.19 update: adding platform_code in dataset information display
% 21.03.23 update: adding different import for platform metadata when
% platform = float
% 21.07.15 update: adding np display in disp(dataset_info_temp)
% 
% -----------------------------------------------------------------------

disp(strcat('loading DATA for platform:',...
    tagRef,'.....'))


%% load dataset

switch platform_type
    case 'sealtag'
	% load structure with all fields
	platform = ncload_struct([root_data tagRef]) ;
    case 'float'
        % load structure with all fields
        platform = ncload_struct([root_data tagRef]) ;
    otherwise
        disp('WARNING: platform is neither sealtag nor float')
end

% display platform dataset properties and platform metadata
ncdisp([root_data tagRef]) ;


%% all datasets do not have LATITUDE_ADJUSTED/LONGITUDE_ADJUSTED fields
% if present, erase LATITUDE/LONGITUDE fields by corresponding *_ADJUSTED fields

if isfield(platform,'LATITUDE_ADJUSTED') == 0
else
    platform.LATITUDE = platform.LATITUDE_ADJUSTED ;
    platform.LONGITUDE = platform.LONGITUDE_ADJUSTED ;
end


%% clear temp data
clear *_temp

%% END