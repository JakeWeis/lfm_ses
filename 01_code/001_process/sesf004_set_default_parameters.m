%% script metadata
%
% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% set default parameters for project
%
%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x
%
%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])
%
%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% J. Weis (jw)
% last modified: 22.03.29
% 
% 20.02.10 update:
% adding delta_for_cstLIGHT parameter for use in dark signal preprocessing
% 20.04.27 update: use of N-dimensional find to determine max depth of dataset
% N-dimensional find
% version 1.2.0.0 (1.99 KB) by Sven
% https://www.mathworks.com/matlabcentral/fileexchange/27755-n-dimensional-find
% 20.05.19 update:
% adding latLim and LonLim for map plots
% adding variable for bathymetry data file name and path
% adding variable for satellite map file name and path
% 20.07.01 update:
% add setting of value for min model bathy threshold
% (default = -2000 m)
% 20.08.03 update: adding readtable to write path for PCA data
% 20.08.26 update: set maxCHLA_depth to 200 dbar (P. Lovell 20.08.24)
% 20.08.28 udpate: modifying deepest nonNaN value computation and
% implementing truncation of dataset to max_depth user-defined value
% (default = 1000 m)
% 20.08.31 updates:
% modifying deepest nonNaN value computation by using TEMP instead of PRES
% (PRES is sometimes non NaN where TEMP is)
% extending lonLim for DDU from 145º to 150º
% 20.09.08 update: updating computation of max_depth_dataset
% 20.11.23 update: adding linux etopo path for linux os
% 21.01.29 update: adjusting database (bathy/sat) paths
% 21.03.24 update: approximating max_depth_dataset with
% floor(max_depth_dataset)
% 21.04.06 update: modifying sat database path and removing
% platform_metadada.satmapfilename
% 21.04.07 update: setting delta_for_cstCHLA = 0.01 
% 21.04.27 updates:
% changing max_depth warning message
% removing min_integ_interval variable (unused)
% 21.07.16 update: adding fda parameters
% 22.02.01 update: changing PAF latlim from [-60 -35] to [-65 -35]
% 22.02.02 update: add MLD min start depth parameter and min profile length
% 22.02.03 update: add topBound_lfm and botBound_lfm
% 22.02.04 update: adding setting of nProfilesPerDay_dark (used in light
% dark signal correction in sesf032)
% 22.02.05 update: adding minProfDepthOpenOcean (default = - 500 m) to
% differentiate open ocean vs coastal when no bathymetry data available
% 22.10.28 updates:
% adding OcCorrFactorMerged = 5.942344745608111 correction
% factor based on satellite regression (all ct152+ct159+ft22 tags merged)
% adding Linear Functional Model (LFM) default parameters
% -----------------------------------------------------------------------

disp(strcat('setting DEFAULT parameters for platform:',...
    platform_metadata.platform_name,'.....'))

%% general default parameters

% total number of profiles
np_tot = length(platform.JULD) ; % total number of profiles in dataset


% max depth to consider for overall dataset
PRES_temp = platform.PRES(:) ;
% Find k largest elements of platform.PRES
[max_temp,ind_temp] = maxk(PRES_temp,length(PRES_temp)) ;
TEMP_temp = platform.TEMP(:) ;
% initialize a few variables
max_depth_dataset = max(PRES_temp) ;
flag_temp = 0 ;
ii_temp = 0 ;
% Find the largest elements of platform.PRES where platform.TEMP is not NaN
while ii_temp < numel(ind_temp)
    ii_temp = ii_temp + 1 ;
    i_temp = ind_temp(ii_temp) ;
    max_depth_dataset = PRES_temp(i_temp) ;
    if isnan(TEMP_temp(i_temp)) == 0
        break
    end
end
% convert max_depth_dataset into integer
max_depth_dataset = floor(max_depth_dataset) ;


if isnan(max_depth_dataset) | max_depth_dataset == 1
    disp('WARNING: looks like dataset is empty')
end

max_depth = 1000 ; % set default to 1000 m
disp(strcat('deepest record in dataset:',{' '},int2str(max_depth_dataset),'m'))
disp(strcat('user-defined deepest value:',{' '},int2str(max_depth),'m'))


% other depth parameters
min_topIntegBound = 50 ;    % minimum upper boundary (surface) required for profiles to be imported in linear functional model
min_botIntegBound = 200 ;   % minimum lower boundary (bottom) required for profiles to be imported in linear functional model

% minimum profile depth
minProfileDepth = 20 ; % (meters)
minProfileLength = 20 ; % (meters)

% DO/DO NOT display plots in preprocessing
doPlots = 'YES' ;

% theta s diagram parameter
tsDiagDepth = 200 ; % depth at which pick up T and S values for T S Diagram


%% default parameters for CHLA

maxCHLA_depth = 200 ;      % define depth at which to consider dark signal (CHLA "absolute zero")
delta_for_cstCHLA = 0.01 ;  % define standard delta between min and max of profile to eliminate profile for being a constant profile


%% default parameters for LIGHT

% use of solar position function
TimOffZ = 0;                 % [hrs] offset from UTC, during standard time
% Date and time values are always in universal time coordinates (UTC)
% see seamammal_user_manual_version1.0.pdf
DaySavTim = false;            % local time is without daylight savings time

% daylight filter
sol_elv_threshold = 0 ;	% threshold to be applied for daylight filter (>20º solar elevation)
delta_for_cstLIGHT = 1E-8 ;  % define standard delta between min and max of profile to eliminate profile for being a constant profile
Epsilon = 1.0000e-6;	% set Epsilon>0 for PAR negative values (further use of log)

% dark correction
% select a few profiles per day for dark correction: the nProfiles deepest per day
nProfilesPerDay_dark = 5 ;


%% default parameters for functional data analysis
data_temp = load([root_data_seal 'LFMcore.mat']) ;
LFMcore = data_temp.LFMcore ;
nbaz = 30 ; % number of basis functions in functional space
nord = 4 ; % order of B-splines
topBound_lfm = 5 ;
botBound_lfm = 200 ;
lambda = 0.03 ; % smoothing parameter for functional fit
mybL = LFMcore.mybL ;
pp = (topBound_lfm:1:botBound_lfm).' ;

%% logical indexing to manually discard profiles along data processing
% (ex. if test/not consistent)

idx000_genProfileSelec = true(np_tot,1) ;


%% logical indexing for non NaN LAT/LON

idx001_latlonNonNan = and(~isnan(platform.LATITUDE),~isnan(platform.LONGITUDE)) ;
% idx001_latlonNonNan = idx001_latlonNonNan.' ; % transpose matrix to get a 1 x n array (display homogeneity)


%% base station location for geo plots
switch base_station
    case 'DDU'
    platform_metadata.base_station_LAT = -66.663253 ;
    platform_metadata.base_station_LON = 140.002335 ;
    platform_metadata.latLim = [-70 -65] ;
    platform_metadata.lonLim = [135 150] ;
    case 'PAF'
    platform_metadata.base_station_LAT = -49.351679 ;
	platform_metadata.base_station_LON = 70.219615 ;
	platform_metadata.latLim = [-65 -35] ;
    platform_metadata.lonLim = [35 100] ;
    case 'P.Del.'
	platform_metadata.base_station_LAT = -42.767604 ;
	platform_metadata.base_station_LON = -63.635241 ;
    platform_metadata.latLim = [-60 -35] ;
    platform_metadata.lonLim = [-80 -40] ;
    otherwise      
    disp('WARNING: base station is nor PAF nor P.Del. nor DDU')
end

%% PCA (prey capture atteimpts) data

% read table containing sealtag naming ocrrespondances
SEAL_REF  = readtable([root_data_seal 'SEAL_REF.csv']) ;
% extract name of platform for PCA data
try
    animalId = SEAL_REF.id_animal(strcmp(SEAL_REF.ref_smru,...
        platform_metadata.platform_name) == true) ;
catch
    animalId = NaN ;
end
platform_metadata.pcaDataPath = root_data_seal ;


%% Bathymetry database to be used

platform_metadata.bathyDataPath =...
    [root_proj '00_data' filesep '00_sat_maps_coastline_bathymetry' filesep 'ETOPO1' filesep] ;
platform_metadata.bathyDataFileName =...
    'etopo1_bed_c_f4.flt' ;
% ETOPO1
% https://www.ngdc.noaa.gov/mgg/global/
% Cite ETOPO1: doi:10.7289/V5C8276M

% set min bathy threshold for Chl-a/light model to be run on profile
% (case 1 waters e.g. open ocean)
bathyMinModelThreshold = -1500 ; % (meters)
% criteria relative to profile depth to consider open ocean vs coastal
% when no bathymetry data available (ex. when no loc)
minProfDepthOpenOcean = - 500 ; % (meters)


%% satellite map to be used

platform_metadata.satMapPath =...
    [root_proj '00_data' filesep '00_sat_maps_coastline_bathymetry' filesep 'CMEMS' filesep] ;

OcCorrFactorMerged = 5.942344745608111 ;


%% clear temp variables
clear *_temp


%% END