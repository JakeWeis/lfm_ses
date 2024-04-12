%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% save workspace

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
% last modified: 22.04.06
% 
% -----------------------------------------------------------------------

disp(strcat('saving WORKSPACE for platform:',...
    platform_metadata.platform_name,'.....'))

%% save workspace

save(strcat(root.data.workspace,...
    datestr(date,'yymmdd'),...
    '_', platform_metadata.platform_name)) ;

%% END

