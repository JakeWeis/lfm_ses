%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% load platform metadata

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define platform.m
% sesf00x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% last modified: 22.03.29
% 
% 20.09.30 updates:
% removing definition of working folders (now in sesf000_define_platform)
% adding quick display of dataset information
% 20.11.19 update: adding platform_code in dataset information display
% 21.03.23 update: adding different import for platform metadata when
% platform = float
% 21.07.15 update: adding np display in disp(dataset_info_temp)
% 22.02.03 update:
% add writing of smru_platform_name from tagRef case does not exist
% idem for wmo (default:smru code)
% idem for platform number
% 
% -----------------------------------------------------------------------

disp(strcat('loading METADATA for platform:',...
    tagRef,'.....'))


%% load dataset

switch platform_type
    case 'sealtag'
        % load tag metadata
        platform_metadata = ncloadatt_struct([root_data tagRef]) ;
    case 'float'
        % load platform metadata
        tagMetaRef_temp = strrep(tagRef,'Sprof','meta') ;
        platform_metadata = ncload_struct([root_data tagMetaRef_temp]) ;
    otherwise
        disp('WARNING: platform is nor sealtag nor float')
end

% display platform dataset properties and platform metadata
ncdisp([root_data tagRef]) ;



%% write platform name/platform code in platform_metadata

% smru platform name
if isfield(platform_metadata,'smru_platform_code')
else
    platform_metadata.smru_platform_code =...
        erase(tagRef,strcat('_',dataset_type,int2str(1),'_prof.nc')) ;
end

% wmo
if isfield(platform_metadata,'wmo_platform_code')
else
    platform_metadata.wmo_platform_code =...
        platform_metadata.smru_platform_code ;
end

% platform code
if isfield(platform_metadata,'platform_code')
else
    platform_metadata.platform_code =...
        platform.PLATFORM_NUMBER ;
end


%% write platform type in platform_metadata

switch platform_type
    case 'sealtag'
        % platform name
        platform_metadata.platform_name = platform_metadata.smru_platform_code ;
        % platform type
        platform_metadata.platform_type = platform_type ;
        % base station
        platform_metadata.base_station = base_station ;
    case 'float'
        % platform name
        platform_name_temp = convertCharsToStrings(platform.PLATFORM_NUMBER(:,1)) ;
        % remove spaces in name
        platform_name_temp = strrep(platform_name_temp,' ','') ;
        % newStr = strrep( str , old , new )
        % replaces all occurrences of old in str with new
        platform_metadata.platform_name = platform_name_temp ;
        % platform type
        platform_metadata.platform_type = platform_type ;
        % base station
        platform_metadata.base_station = base_station ;
    otherwise
        disp('WARNING: platform is nor sealtag nor float')
end


%% number of profiles / number of rows
platform_metadata.np = length(platform.JULD) ;
platform_metadata.nr = length(platform.TEMP(:,1)) ;


%% WMO / SMRU Nº
switch platform_type
    case 'float'
        wmoNames_temp = horzcat(platform_metadata(1:end).PLATFORM_NUMBER) ;
    case 'sealtag'
        wmoNames_temp = vertcat(platform_metadata(1:end).wmo_platform_code).' ;
        % SMRU = SEALTAG ONLY
        smruNames_temp = vertcat(platform_metadata(1:end).smru_platform_code).' ;
        % create char array of proper size for platform ref
        smru_platform_code_temp = char(platform.TEMP(1:length(smruNames_temp(:,1)),:)) ;
        % write proper ref number for each profile
        smru_platform_code_temp = repmat(smruNames_temp.',platform_metadata.np,1).' ;

        platform_metadata.smru_platform_code = smru_platform_code_temp ;
        platform_metadata.smru_list = cellstr(smruNames_temp.') ;
end

% WMO is both for FLOAT + SEALTAG
% create char array of proper size for ref
wmo_platform_code_temp = char(platform.TEMP(1:length(wmoNames_temp(:,1)),:)) ;
% write proper ref number for each profile
wmo_platform_code_temp = repmat(wmoNames_temp.',platform_metadata.np,1).' ;
platform_metadata.wmo_platform_code = wmo_platform_code_temp ;
platform_metadata.wmo_list = cellstr(wmoNames_temp.') ;

%% write comment in platform_metadata to differentiate weather
% merged/single platform
platform_metadata.comment = 'single_platform' ;


%% additional metadata
    % number of profiles
platform_metadata.number_of_ts_profiles = platform_metadata.np ;
platform_metadata.number_of_t_profiles = platform_metadata.np ;
platform_metadata.number_chla_profiles = platform_metadata.np ;
platform_metadata.number_light_profiles = platform_metadata.np ;
    % lat/lon min/max
platform_metadata.geospatial_lat_max = max(platform.LATITUDE) ;
platform_metadata.geospatial_lat_min = min(platform.LATITUDE) ;
platform_metadata.geospatial_lon_max = max(platform.LONGITUDE) ;
platform_metadata.geospatial_lon_min = min(platform.LONGITUDE) ;


%% diplay information
dataset_info_temp = { platform_metadata.platform_name...
    base_station platform_type dataset_type length(platform.JULD)} ;
dataset_info_temp = cell2table(dataset_info_temp,...
    'VariableNames', {'platform_code' 'base_station' 'platform_type' 'dataset_type' 'n profiles'}) ;
disp(dataset_info_temp)


%% now ready to execute preprocessing steps...


%% clear temp data
clear *_temp

%% END