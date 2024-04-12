%% lfm_ses_menu
%
% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% This set of routines predicts Chla from Kd from marine mammal data
% This routines controls everything, it
% imports the tag data,
% runs data pre-processing,
% runs Linear Functional Model prediction
% saves the output
%
%%%%%%%%%% REFERENCES %%%%%%%%%%
%
%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% none
%
%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])
%
%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% J. Weis (jw)
% last modified: 22/03/2024
% 
% 22/03/2024 update: 
% - Outsourced input data (file sizes too large to push to GitHub)
% - Removed ETOPO file and disabled bathymetry preprocessing steps for the same reason
% -------------------------------------------------------------------------

%% clear all and define project root/input/output directories
% Get root directory and add to MATLAB path
root.proj = mfilename('fullpath');
i_filesep = strfind(root.proj,filesep);
root.proj(i_filesep(end-1)+1:end) = [];
addpath(genpath(root.proj))

root.data.seal      = [root.proj '00_data' filesep '01_seal' filesep];
root.data.float     = [root.proj '00_data' filesep '02_float' filesep];
root.data.cruise    = [root.proj '00_data' filesep '03_cruise' filesep];
root.data.workspace = [root.proj '00_data' filesep '04_workspace' filesep];
root.plots          = [root.proj '06_plots' filesep];

% input data directory (TO BE SPECIFIED AS INPUT TO)
% root.input = '/Volumes/PhData/PD DATA/MEOP-CTD_2024-03-08/SUBSET/FLUO_LIGHT';
root.input          = '/Volumes/PhData/PD DATA/test_data/';
if ~strcmp(root.input(end),filesep)
    root.input      = [root.input filesep];
end
root.output         = [root.input 'LFM_SES_output' filesep];
if ~isfolder(root.output)
    mkdir(root.output)
end


%% display tag names

sesf000_display_tags_names
% F000_DisplayTagNames(root)

%% Load ETOPO Global Relief Model
% NOAA National Centers for Environmental Information (2022). ETOPO 2022 60 Arc-Second Global Relief Model.
% DOI: https://doi.org/10.25921/fd45-gt74. Accessed 19/03/2024.
if ~exist('bathymetry','var')
    disp('Loading ETOPO 2022 Global Relief Model...')
    [bathymetry.data, bathymetry.ref] = readgeoraster([root.proj, '00_data', filesep, '00_etopo', filesep, 'ETOPO_2022_v1_60s_N90W180_bed.tif']);
    bathymetry.data = flipud(bathymetry.data);
    bathymetry.lon = bathymetry.ref.LongitudeLimits(1) : bathymetry.ref.CellExtentInLongitude : bathymetry.ref.LongitudeLimits(2) - bathymetry.ref.CellExtentInLongitude;
    bathymetry.lat = bathymetry.ref.LatitudeLimits(1) : bathymetry.ref.CellExtentInLatitude : bathymetry.ref.LatitudeLimits(2) - bathymetry.ref.CellExtentInLatitude;
end

%% data processing steps

for iTag = 1 : numel(fold_info)

    
	%% clean up workspace
    close all
    clearvars -except iTag fold_info* root*...
        platform_type base_station dataset_type...
        sep bathymetry
    
    %% load platform
    
    tagRef = fold_info(iTag).name ;
	disp(char(strcat('tag:',{' '},tagRef)))
	% load corresponding dataset
    sesf003_load_data_single_platform
    sesf003_load_metadata_single_platform    
    
    
    %% APPETIZERS: set default parameters

    % set default parameters
    sesf004_set_default_parameters

    % reedit parameters
    doPlots = 'NO' ; % default: doPlots = 'YES' ;


    %% STARTERS: date, pressure, trip distance, set regular grid and explore bathymetry
    
    sesf010_set_date_arrays_genprofdata
    
    % eventually remove some profiles (if correspond to test)
    % e.g. remove test profiles performed in August (month = 8)
    monthsToRemove_temp = [2 3 4 5 6 7 8] ;
    toBeRemoved_temp = arrayfun(@(a) month(dateYMD(1,:)) == a, monthsToRemove_temp,'UniformOutput',false) ;
    toBeRemoved_temp = vertcat(toBeRemoved_temp{1:end}) ;
    toBeRemoved_temp = sum(toBeRemoved_temp) ;
    toBeRemoved_temp = logical(toBeRemoved_temp) ;
    toBeRemoved_temp = [] ;
    idx000_genProfileSelec(toBeRemoved_temp) = false ;
    genData.depDayNo(toBeRemoved_temp) = NaN ;
    genData.depDayNo = genData.depDayNo - min(genData.depDayNo) + 1 ;

    
    sesf011_gen_data_preprocess1_trip_distance
    sesf012_gen_data_preprocess2_set_fixed_vert_grid
    sesf013_gen_data_preprocess3_bathy


    %% WITH A DRINK: a few talks about density

    sesf015_gen_data_preprocess5_compute_dens
    sesf018_gen_data_preprocess7_MLDphy003
        genData.MLDHTd2009 = nan(np_tot,1) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % choose which MLD computation you want to apply (if multiple methods
    % were used)
    % (MLD003/MLDphyBw/MLDphyDHT2009/MLDphyTHT2009/MLDPhySlm/MLDt02)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    genData.MLDphy = fillmissing(genData.MLD003,'nearest','EndValues','none') ;


    %% MAIN COURSE: light

    sesf031_par_data_preprocess1_nonNan
    sesf032_par_data_preprocess2_dark_signal
    sesf033_par_data_preprocess3_saturation
    sesf036_par_data_preprocess6_functionalFit
    sesf037_par_data_preprocess7_postMetrics

    %% SALAD: fluorescence

    sesf040_fluo_data_preprocess0_nonNan
    sesf041_fluo_data_preprocess1_rawMetrics
    sesf042_fluo_data_preprocess2_darkSignal
    sesf043_fluo_data_preprocess3_RSD_mixingLayer
    sesf044_fluo_data_preprocess4_NPQX18
    sesf046_fluo_data_preprocess6_smoothFluo
    sesf047_fluo_data_preprocess7_postMetrics


    %% LITTLE BREAK: admire plots


    %% TIME FOR DESERT: create logical indexes for custom export
    % general logical index (LAT/LON/CHLA/LIGHTH non NaN)
    sesf100_create_logical_index
    % custom export logical index
    sesf101_create_custom_logical_index

    
    %% DESERT: check out satellite data
    
    
    %% FRUIT: chlorophyll-a
 
    % predict Chl-a from Kd profiles
    sesf270_chla_LFMprediction
    
    % chlaLFM post-processing metrics
    sesf347_chlaLFM_data_postMetrics
       
        
    %% still some room for additional analysis?
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % customize here your analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % end of custom analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% THE BILL: save output, workspace and plots
    sesf997_save_output
    % sesf998_save_workspace
%     sesf999_save_plots

end



%% clear temp variables
clear *_temp


%% END