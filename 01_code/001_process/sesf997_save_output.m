%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% save LFM prediction

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% J. Weis (jw)
% last modified: 22.10.31
% 
% -----------------------------------------------------------------------

disp(strcat('saving OUTPUT for platform:',...
    platform_metadata.platform_name,'.....'))

%% save workspace

save(strcat(root_data_output,...
    datestr(date,'yymmdd'),...
    '_CHLALFM_', platform_metadata.platform_name,'.mat'),'CHLA_LFM') ;


%% END

