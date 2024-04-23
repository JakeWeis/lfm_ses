function [platform_data,platform_metadata,platform_processed,genData] = loadData(root,tagRef,defVals,bathymetry)
% LOAD DATA loads the seal tag data and metadata.
%
% INPUT ARGUMENTS
%
% OUTPUT

%% CMD message: start
fprintf('Loading data...');

%% Identify data type from NetCDF attributes (animal-borne or Argo)
% Might remove the Argo bit at some point as this code is meant to be seal
% tag specific...
if strcmp(ncreadatt([root.input tagRef],'/','source'), 'Marine mammal observation')
    platform_type = 'sealtag';
elseif strcmp(ncreadatt([root.input tagRef],'/','source'), 'Argo float')
    platform_type = 'float';
else
    platform_type = 'undefined';
    warning('File is neither seal tag nor Argo float.')
end

%% Load data
platform_data = ncload_struct([root.input tagRef]);

% If available, replace LATITUDE/LONGITUDE fields with corresponding *_ADJUSTED fields
if isfield(platform_data,'LATITUDE_ADJUSTED')
    platform_data.LATITUDE = platform_data.LATITUDE_ADJUSTED;
    platform_data.LONGITUDE = platform_data.LONGITUDE_ADJUSTED;
end

% Convert JULD to DATETIME and DATENUM
reference_date = '1950-01-01 00:00:00'; % see platform_data.REFERENCE_DATE_TIME
dateYMD = datetime(platform_data.JULD_LOCATION + datenum(reference_date),'ConvertFrom','datenum');
platform_data.DATETIME = repmat(dateYMD',size(platform_data.PRES,1),1);

% Calculate depth from pressure and latitude (using GSW equations)
platform_data.DEPTH = gsw_z_from_p(platform_data.PRES_ADJUSTED,platform_data.LATITUDE);

% Calculate density from abs. salinity, conservative temperature and pressure (using GSW equations)
SA = gsw_SA_from_SP(platform_data.PSAL_ADJUSTED,platform_data.PRES_ADJUSTED,platform_data.LONGITUDE,platform_data.LATITUDE);
CT = gsw_CT_from_t(SA,platform_data.TEMP_ADJUSTED,platform_data.PRES_ADJUSTED);
platform_data.DENS = gsw_rho_CT_exact(SA,CT,0);

%% Load metadata
switch platform_type
    case 'sealtag'
        % Seal tag metadata stored in main data file
        platform_metadata = ncloadatt_struct([root.input tagRef]);
        platform_metadata.platform_type = 'sealtag';
    case 'float'
        % Float metadata stored in separate file
        float_metadata_filename = strrep(tagRef,'Sprof','meta');
        platform_metadata = ncload_struct([root.input float_metadata_filename]);
        platform_metadata.platform_type = 'float';
end

% Number of profiles and number of depth levels
platform_metadata.nProfs = size(platform_data.DEPTH,2);
platform_metadata.nDepth = size(platform_data.DEPTH,1);

% Check resolution of tag data
if ~isempty(regexp(tagRef,'lr[1-9]','once')) || platform_metadata.nDepth < numel(defVals.depthInterpGrid)
    % if file name contains "lr1/lr2/..."  or 
    % if n depth levels < n depth levels of the specified depth interpolation vector
    platform_metadata.depth_res = 'lr';
elseif ~isempty(regexp(tagRef,'hr[1-9]','once')) || platform_metadata.nDepth == numel(defVals.depthInterpGrid)
    % if file name contains "hr1/hr2/..." or i
    % if n depth levels = n depth levels of the specified depth interpolation vector
    platform_metadata.depth_res = 'hr';
elseif ~isempty(regexp(tagRef,'fr[1-9]','once'))
    % if file name contains "fr1/fr2/..."
    platform_metadata.depth_res = 'fr';
end

% lat/lon min/max
platform_metadata.geospatial_lat_max = max(platform_data.LATITUDE);
platform_metadata.geospatial_lat_min = min(platform_data.LATITUDE);
platform_metadata.geospatial_lon_max = max(platform_data.LONGITUDE);
platform_metadata.geospatial_lon_min = min(platform_data.LONGITUDE);

% sort depth from deepest to shallowest and find the deepest depth where temperature is not NaN
PRES_vec = platform_data.PRES(:);
[~,i_sort] = sort(PRES_vec,'descend','MissingPlacement','last');
TEMP_vec = platform_data.TEMP(:);
platform_metadata.max_depth_dataset = max(PRES_vec);
for iP = 1 : numel(i_sort)
    platform_metadata.max_depth_dataset = PRES_vec(i_sort(iP));
    if isfinite(TEMP_vec(i_sort(iP)))
        break
    end
end

%% Interpolate data onto uniform depth grid

switch platform_metadata.depth_res

    case 'lr' % low resolution

        % The depth array of the processed data is equivalent to the predefined interpolation-depth vector (in meters)
        platform_processed.DEPTH = repmat(defVals.depthInterpGrid,1,platform_metadata.nProfs);
        
        % Pressure, temperature, salinity, fluorescence and PAR are all interpolated onto the uniform depth grid given above
        platform_processed.PRES.Reg = NaN(numel(defVals.depthInterpGrid),platform_metadata.nProfs);
        platform_processed.TEMP.Reg = NaN(numel(defVals.depthInterpGrid),platform_metadata.nProfs);
        platform_processed.PSAL.Reg = NaN(numel(defVals.depthInterpGrid),platform_metadata.nProfs);
        platform_processed.FLUO.Reg = NaN(numel(defVals.depthInterpGrid),platform_metadata.nProfs);
        platform_processed.PAR.log_Reg = NaN(numel(defVals.depthInterpGrid),platform_metadata.nProfs);
        platform_processed.PAR.lin_Reg = NaN(numel(defVals.depthInterpGrid),platform_metadata.nProfs);

        % Interpolate
        for iP = 1 : platform_metadata.nProfs
            % depth vector for profile i (irregular step)
            depth_irreg = platform_data.DEPTH(:,iP);
            depth_isnan = ~isfinite(depth_irreg);
            depth_irreg(depth_isnan) = [];

            % Pressure
            pres_irreg = platform_data.PRES_ADJUSTED(:,iP);
            pres_irreg(depth_isnan) = [];
            platform_processed.PRES.Reg(:,iP) = interp1(depth_irreg,pres_irreg,defVals.depthInterpGrid,'linear');

            % Temperature
            temp_irreg = platform_data.TEMP_ADJUSTED(:,iP);
            temp_irreg(depth_isnan) = [];
            platform_processed.TEMP.Reg(:,iP) = interp1(depth_irreg,temp_irreg,defVals.depthInterpGrid,'linear');

            % Practical salinity
            psal_irreg = platform_data.PSAL_ADJUSTED(:,iP);
            psal_irreg(depth_isnan) = [];
            platform_processed.PSAL.Reg(:,iP) = interp1(depth_irreg,psal_irreg,defVals.depthInterpGrid,'linear');

            % Chl fluorescence
            if isfield(platform_data,'CHLA')
                fluo_irreg = platform_data.CHLA(:,iP);
                fluo_irreg(depth_isnan) = [];
                platform_processed.FLUO.Reg(:,iP) = interp1(depth_irreg,fluo_irreg,defVals.depthInterpGrid,'linear');
            end

            % LIGHT fluorescence
            if isfield(platform_data,'LIGHT')
                par_irreg = platform_data.LIGHT(:,iP);
                par_irreg(depth_isnan) = [];
                platform_processed.PAR.log_Reg(:,iP) = interp1(depth_irreg,par_irreg,defVals.depthInterpGrid,'linear');
                platform_processed.PAR.lin_Reg(:,iP) = exp(platform_processed.PAR.log_Reg(:,iP));
            end
        end

    otherwise % high or full resolution

        platform_processed.DEPTH = platform_data.DEPTH;
        platform_processed.PRES.Reg = platform_data.PRES_ADJUSTED;
        platform_processed.TEMP.Reg = platform_data.TEMP_ADJUSTED;
        platform_processed.PSAL.Reg = platform_data.PSAL_ADJUSTED;
        platform_processed.FLUO.Reg = platform_data.CHLA;
        platform_processed.PAR.log_Reg = platform_data.LIGHT;
        platform_processed.PAR.lin_Reg = exp(platform_data.LIGHT);

end

%% create genData table to store general profile information
var_names = {...
    'Profile', ...                  % profile number
    'TagID', ...                    % seal tag ID
    'Date', ...                     % date
    'DeployDay', ...                % deployment day
    'Lon', ...                      % lon
    'Lat', ...                      % lat
    'Distance', ...                 % distance travelled by seal between profiles
    'CumDistance', ...              % cumulative distance travelled by seal across all profiles
    'Processed', ...                % true/false: profile processed
    'BathymetryAtProf', ...         % bathymetry at profile location
    'TimeOfDay_UTC', ...            % hour of the day (UTC) of the profile
    'SolarAlt', ...                 % solar altitude value (angle, in dregrees)
    'Daytime', ...                  % true/false: daytime profile
    'MLD', ...                      % MLD [T, S, RHO] (based on Holte and Talley 2009 algorithms)
    }';

genData = array2table(NaN(platform_metadata.nProfs, numel(var_names)),'VariableNames',var_names);

% Write data to genData table
genData.Profile = (1 : platform_metadata.nProfs)';
genData.TagID = repmat(platform_metadata.platform_code,platform_metadata.nProfs,1);
genData.Date = dateYMD;
genData.DeployDay = floor(datenum(dateYMD)) - floor(datenum(dateYMD(1))) + 1;
genData.Lon = platform_data.LONGITUDE;
genData.Lat = platform_data.LATITUDE;

% Shallowest and deepest available observation (based on temperature)
firstNonNaN = NaN(platform_metadata.nProfs,2);
lastNonNaN = NaN(platform_metadata.nProfs,2);
firstNonNaN(:,1) = find_ndim(isfinite(platform_data.TEMP_ADJUSTED),1,'first').';
lastNonNaN(:,1) = find_ndim(isfinite(platform_data.TEMP_ADJUSTED),1,'last').';
for iP = 1 : platform_metadata.nProfs
    if firstNonNaN(iP,1) > 0
        firstNonNaN(iP,2) = platform_data.DEPTH(firstNonNaN(iP,1),iP);
    end
    if lastNonNaN(iP,1) > 0
        lastNonNaN(iP,2) = platform_data.DEPTH(lastNonNaN(iP,1),iP);
    end
end

% Exclude profiles from processing where no usable observations are available
genData.Processed = true(platform_metadata.nProfs,1);
genData.Processed(firstNonNaN(:,2) < defVals.depthInterpGrid(end)) = false; % no obs above maximum profile depth
genData.Processed(lastNonNaN(:,2) >= defVals.minProfileDepth) = false; % no obs below minimum default profile depth
genData.Processed(firstNonNaN(:,2) - lastNonNaN(:,2) < defVals.minProfileLength) = false; % depth range of usable obs < minimum default profile length

% Distance travelled by seal
distance_km = zeros(platform_metadata.nProfs,2);
for iP = 2 : platform_metadata.nProfs
    if all(isfinite([platform_data.LATITUDE(iP-1:iP);platform_data.LONGITUDE(iP-1:iP)]))
        loc1 = [platform_data.LATITUDE(iP-1) platform_data.LONGITUDE(iP-1)] ;
        loc2 = [platform_data.LATITUDE(iP) platform_data.LONGITUDE(iP)] ;
        [dist1km_temp dist2km_temp] = lldistkm(loc1, loc2) ;
        % write computed values in platform_distancekm array
        distance_km(iP,1) = distance_km(iP-1,1) + dist1km_temp ;  % col 1 = cumulated distance
        distance_km(iP,2) = dist1km_temp ;	% col 2 = relative distance
    else
        distance_km(iP,1) = distance_km(iP-1,1) ;
        distance_km(iP,2) = NaN ;
    end
end
genData.Distance = distance_km(:,2);
genData.CumDistance = distance_km(:,1);

%% Solar altitude
dateNUM = datenum(dateYMD);
for iP = 1 : platform_metadata.nProfs
        % compute solar angle and solar altitude (even for all NaN profiles)   
        % solar angle
        zenith_azimuth = solarPosition( ...
            dateNUM(iP),...
            platform_data.LATITUDE(iP), ...
            platform_data.LONGITUDE(iP),...
            defVals.PAR.UTC_offset,0,defVals.PAR.DST);

    	% write computed values in genData array
        genData.SolarAlt(iP) = 90 - zenith_azimuth(1); % solar altitude value (angle, in dregrees)
    	genData.TimeOfDay_UTC(iP) = datepart(dateYMD(iP),'hour'); % hour of the day (UTC) of the profile
end

% Identify daytime profiles
genData.Daytime = true(platform_metadata.nProfs,1);
genData.Daytime(genData.SolarAlt <= defVals.PAR.SolarAlt_DayTime) = false;

%% Get bathymetry data along seal path
% Extract bathymetry subset around the seal track's lat/lon limits from ETOPO dataset
i_latmin_in_bathy = dsearchn(bathymetry.lat',floor(platform_metadata.geospatial_lat_min));
i_latmax_in_bathy = dsearchn(bathymetry.lat',ceil(platform_metadata.geospatial_lat_max));
i_lonmin_in_bathy = dsearchn(bathymetry.lon',floor(convlon(platform_metadata.geospatial_lon_min,'signed')));
i_lonmax_in_bathy = dsearchn(bathymetry.lon',ceil(convlon(platform_metadata.geospatial_lon_max,'signed')));

if platform_metadata.geospatial_lon_max - platform_metadata.geospatial_lon_min < 180
    bathy_subset = bathymetry.data(i_latmin_in_bathy:i_latmax_in_bathy, i_lonmin_in_bathy:i_lonmax_in_bathy);
    bathy_lon = bathymetry.lon(i_lonmin_in_bathy : i_lonmax_in_bathy);
    bathy_lat = bathymetry.lat(i_latmin_in_bathy : i_latmax_in_bathy);
else
    % If the lon limits cross 180˚ E/W, the bathymetry subset needs to be pieced together
    bathy_E = bathymetry.data(i_latmin_in_bathy:i_latmax_in_bathy, i_lonmax_in_bathy:bathymetry.ref.RasterSize(2));
    bathy_W = bathymetry.data(i_latmin_in_bathy:i_latmax_in_bathy, 1:i_lonmin_in_bathy);
    bathy_subset = [bathy_E,bathy_W];
    bathy_lon_E = bathymetry.lon(i_lonmax_in_bathy:bathymetry.ref.RasterSize(2));
    bathy_lon_W = bathymetry.data(1:i_lonmin_in_bathy);
    bathy_lon = [bathy_lon_E,bathy_lon_W];
    bathy_lat = bathymetry.lat(i_latmin_in_bathy : i_latmax_in_bathy);
end

% Get bathymetry depth at each profile
% Find corresponding seal lat/lon indices in bathymetry subgrid
i_lat_in_bathy = dsearchn(bathy_lat',platform_data.LATITUDE);
i_lat_in_bathy(~isfinite(platform_data.LATITUDE)) = NaN;
i_lon_in_bathy = dsearchn(bathy_lon',platform_data.LONGITUDE);
i_lon_in_bathy(~isfinite(platform_data.LONGITUDE)) = NaN;
genData.BathymetryAtProf = NaN(platform_metadata.nProfs,1);
for iP = 1 : platform_metadata.nProfs
    if isfinite(i_lat_in_bathy(iP)) && isfinite(i_lon_in_bathy(iP))
        genData.BathymetryAtProf(iP) = bathy_subset(i_lat_in_bathy(iP),i_lon_in_bathy(iP));
    end
end

%% MLD
MLD_algo = 3; % 1: Temperature threshold, 2: Salinity threshold, 3: Density threshold
for iP = 1 : platform_metadata.nProfs
    pres = platform_processed.PRES.Reg(:,iP);
    temp = platform_processed.TEMP.Reg(:,iP);
    psal = platform_processed.PSAL.Reg(:,iP);
    mld_estimates_p = mld(pres(isfinite(pres)),temp(isfinite(pres)),psal(isfinite(pres)));
    
    mld_estimates = gsw_z_from_p(mld_estimates_p,repmat(genData.Lat(iP),1,3));

    genData.MLD(iP) = mld_estimates(:,MLD_algo);
end

%% CMD message: done
fprintf('\b\b \x2713\n')

end