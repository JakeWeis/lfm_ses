%% script metadata
% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% compute sea water density using TEMP and PSAL
% functions
%   sw_pden from SeaWater library of EOS-80
%   source: http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf001_choose_your_menu
% sesf003_load_dataset
% sesf004_set_default_parameters
% sesf010_set_date_pressure_arrays
% sesf012_gen_data_preprocess2_set_fixed_vert_grid

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% last modified: 22.04.18
% 
% 20.09.25 update: removing for loop
% 22.04.18 update: adding computation of DENS_RAW (NON-ADJUSTED)
%
% -----------------------------------------------------------------------

disp(strcat('computing sea water DENSITY for platform:',...
    platform_metadata.platform_name,'.....'))

%% compute sea water potential density

ref_pres_temp = 0 ; % Reference pressure [db]
DENS_RAW = sw_pden(platform.PSAL(:,:),platform.TEMP(:,:),platform.PRES(:,:),ref_pres_temp) - 1000 ;
DENS = sw_pden(SAL(:,:),TEMP(:,:),PRES(:,:),ref_pres_temp) - 1000 ;


%% CLEAR TEMP VARIABLES
clear *_temp


%% END
