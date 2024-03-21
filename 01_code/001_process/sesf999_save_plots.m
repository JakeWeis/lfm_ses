%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% save plots (all open figures)

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% last modified: 22.04.06
%
% 21.04.27 update: modifying script to print all open figures
% 22.02.03 update: modifying format_temp from -dpng to -djpeg
%
% -----------------------------------------------------------------------

disp(strcat('saving FIGURES for platform:',...
    platform_metadata.platform_name,'.....'))


%% save plots

def_temp = '-r600' ;
format_temp = '-djpeg' ;

figList_temp = findobj(0, 'type', 'figure');

for iFig_temp = 1:numel(figList_temp)
  figHandle_temp = figList_temp(iFig_temp);
  figNum_temp   = num2str(get(figHandle_temp, 'Number'));
  figFileName_temp = strcat(root_plots,...
        datestr(date,'yymmdd'),...
        '_',platform_metadata.platform_name,...
        '_plot_',figNum_temp) ;
  print(strcat('-f',figNum_temp),figFileName_temp,def_temp,format_temp)
end


%% clear temp variables
clear *_temp

%% END