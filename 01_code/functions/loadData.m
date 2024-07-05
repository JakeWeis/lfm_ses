function [Data,ProfileInfo] = loadData(root,platformID,bathymetry,defaultPars)
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
if strcmp(ncreadatt(fullfile(root.input, platformID),'/','source'), 'Marine mammal observation')
    platform_type = 'sealtag';
elseif strcmp(ncreadatt(fullfile(root.input, platformID),'/','source'), 'Argo float')
    platform_type = 'float';
else
    platform_type = 'undefined';
    warning('File is neither seal tag nor Argo float.')
end

%% Load data
Data.Raw = ncload_struct(fullfile(root.input, platformID));

% If available, replace LATITUDE/LONGITUDE fields with corresponding *_ADJUSTED fields
if isfield(Data.Raw,'LATITUDE_ADJUSTED')
    Data.Raw.LATITUDE = Data.Raw.LATITUDE_ADJUSTED;
    Data.Raw.LONGITUDE = Data.Raw.LONGITUDE_ADJUSTED;
end

% Convert JULD to DATETIME and DATENUM
reference_date = '1950-01-01 00:00:00'; % see platform_data.REFERENCE_DATE_TIME
dateYMD = datetime(Data.Raw.JULD_LOCATION + datenum(reference_date),'ConvertFrom','datenum');
Data.Raw.DATETIME = repmat(dateYMD',size(Data.Raw.PRES,1),1);

% Calculate density from abs. salinity, conservative temperature and pressure (using GSW equations)
SA = gsw_SA_from_SP(Data.Raw.PSAL_ADJUSTED,Data.Raw.PRES_ADJUSTED,Data.Raw.LONGITUDE,Data.Raw.LATITUDE);
CT = gsw_CT_from_t(SA,Data.Raw.TEMP_ADJUSTED,Data.Raw.PRES_ADJUSTED);
Data.Raw.DENS = gsw_rho_CT_exact(SA,CT,0);

%% Load metadata
switch platform_type
    case 'sealtag'
        % Seal tag metadata stored in main data file
        Data.MetaData = ncloadatt_struct(fullfile(root.input, platformID));
        Data.MetaData.platform_type = 'sealtag';
    case 'float'
        % Float metadata stored in separate file
        float_metadata_filename = strrep(platformID,'Sprof','meta');
        Data.MetaData = ncload_struct(fullfile(root.input, platformID));
        Data.MetaData.platform_type = 'float';
end

% Number of profiles and number of depth levels
Data.MetaData.nProfs = size(Data.Raw.PRES,2);
Data.MetaData.nDepth = size(Data.Raw.PRES,1);

% Check resolution of tag data
if ~isempty(regexp(platformID,'lr[1-9]','once'))
    % Low res: If the file name contains "lr1/lr2/..."
    Data.MetaData.depth_res = 'lr';
elseif ~isempty(regexp(platformID,'hr[1-9]','once'))
    % High res: If the file name contains "hr1/hr2/..."
    Data.MetaData.depth_res = 'hr';
elseif ~isempty(regexp(platformID,'fr[1-9]','once'))
    % Full Res: If the file name contains "fr1/fr2/..."
    Data.MetaData.depth_res = 'fr';
else
    % If the resolution is not specified in the file name.
    if Data.MetaData.nDepth < numel(defaultPars.depthInterpGrid)
        % Low res: if n depth levels < n depth levels of the specified depth interpolation vector
        Data.MetaData.depth_res = 'lr';
    elseif Data.MetaData.nDepth == numel(defaultPars.depthInterpGrid)
        % if n depth levels = n depth levels of the specified depth interpolation vector
        Data.MetaData.depth_res = 'hr';
    end
end

% lat/lon min/max
Data.MetaData.geospatial_lat_max = max(Data.Raw.LATITUDE);
Data.MetaData.geospatial_lat_min = min(Data.Raw.LATITUDE);
Data.MetaData.geospatial_lon_max = max(Data.Raw.LONGITUDE);
Data.MetaData.geospatial_lon_min = min(Data.Raw.LONGITUDE);

% sort depth from deepest to shallowest and find the deepest depth where temperature is not NaN
PRES_vec = Data.Raw.PRES(:);
[~,i_sort] = sort(PRES_vec,'descend','MissingPlacement','last');
TEMP_vec = Data.Raw.TEMP(:);
Data.MetaData.max_depth_dataset = max(PRES_vec);
for iP = 1 : numel(i_sort)
    Data.MetaData.max_depth_dataset = PRES_vec(i_sort(iP));
    if isfinite(TEMP_vec(i_sort(iP)))
        break
    end
end

%% Interpolate data onto uniform depth grid

% Calculate depth from pressure and latitude (using GSW equations)
Data.Raw.DEPTH = gsw_z_from_p(Data.Raw.PRES_ADJUSTED,Data.Raw.LATITUDE);

% The depth array of the processed data is equivalent to the predefined interpolation-depth vector (in meters)
Data.Processed.DEPTH = repmat(defaultPars.depthInterpGrid,1,Data.MetaData.nProfs);

% Pressure, temperature, salinity, fluorescence and PAR are all interpolated onto the uniform depth grid given above
Data.Processed.PRES.Reg = NaN(numel(defaultPars.depthInterpGrid),Data.MetaData.nProfs);
Data.Processed.TEMP.Reg = NaN(numel(defaultPars.depthInterpGrid),Data.MetaData.nProfs);
Data.Processed.PSAL.Reg = NaN(numel(defaultPars.depthInterpGrid),Data.MetaData.nProfs);
Data.Processed.FLUO.Reg = NaN(numel(defaultPars.depthInterpGrid),Data.MetaData.nProfs);
Data.Processed.PAR.log.Reg = NaN(numel(defaultPars.depthInterpGrid),Data.MetaData.nProfs);
Data.Processed.PAR.lin.Reg = NaN(numel(defaultPars.depthInterpGrid),Data.MetaData.nProfs);

% Interpolate
for iP = 1 : Data.MetaData.nProfs
    % depth vector for profile i (irregular step)
    depth_irreg = Data.Raw.DEPTH(:,iP);
    depth_isnan = ~isfinite(depth_irreg);
    depth_irreg(depth_isnan) = [];

    % Pressure
    pres_irreg = Data.Raw.PRES_ADJUSTED(:,iP);
    pres_irreg(depth_isnan) = [];
    Data.Processed.PRES.Reg(:,iP) = interp1(depth_irreg,pres_irreg,defaultPars.depthInterpGrid,'linear');

    % Temperature
    temp_irreg = Data.Raw.TEMP_ADJUSTED(:,iP);
    temp_irreg(depth_isnan) = [];
    Data.Processed.TEMP.Reg(:,iP) = interp1(depth_irreg,temp_irreg,defaultPars.depthInterpGrid,'linear');

    % Practical salinity
    psal_irreg = Data.Raw.PSAL_ADJUSTED(:,iP);
    psal_irreg(depth_isnan) = [];
    Data.Processed.PSAL.Reg(:,iP) = interp1(depth_irreg,psal_irreg,defaultPars.depthInterpGrid,'linear');

    % Chl fluorescence
    if isfield(Data.Raw,'CHLA')
        fluo_irreg = Data.Raw.CHLA(:,iP);
        fluo_irreg(depth_isnan) = [];
        Data.Processed.FLUO.Reg(:,iP) = interp1(depth_irreg,fluo_irreg,defaultPars.depthInterpGrid,'linear');
    end

    % LIGHT fluorescence
    if isfield(Data.Raw,'LIGHT')
        par_irreg = Data.Raw.LIGHT(:,iP);
        par_irreg(depth_isnan) = [];
        Data.Processed.PAR.log.Reg(:,iP) = interp1(depth_irreg,par_irreg,defaultPars.depthInterpGrid,'linear');
        Data.Processed.PAR.lin.Reg(:,iP) = exp(Data.Processed.PAR.log.Reg(:,iP));
    end
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
    'SeaIceConcAtProf', ...         % sea ice concentrations at profile location
    'TimeOfDay_UTC', ...            % hour of the day (UTC) of the profile
    'SolarAlt', ...                 % solar altitude value (angle, in dregrees)
    'Daytime', ...                  % true/false: daytime profile
    'MLD', ...                      % MLD [T, S, RHO] (based on Holte and Talley 2009 algorithms)
    }';

ProfileInfo.General = array2table(NaN(Data.MetaData.nProfs, numel(var_names)),'VariableNames',var_names);

% Write data to genData table
ProfileInfo.General.Profile = (1 : Data.MetaData.nProfs)';
ProfileInfo.General.TagID = repmat(Data.MetaData.smru_platform_code,Data.MetaData.nProfs,1);
ProfileInfo.General.Date = dateYMD;
ProfileInfo.General.DeployDay = floor(datenum(dateYMD)) - floor(datenum(dateYMD(1))) + 1;
ProfileInfo.General.Lon = Data.Raw.LONGITUDE;
ProfileInfo.General.Lat = Data.Raw.LATITUDE;

% Shallowest and deepest available observation (based on temperature)
firstNonNaN = NaN(Data.MetaData.nProfs,2);
lastNonNaN = NaN(Data.MetaData.nProfs,2);
firstNonNaN(:,1) = find_ndim(isfinite(Data.Raw.TEMP_ADJUSTED),1,'first').';
lastNonNaN(:,1) = find_ndim(isfinite(Data.Raw.TEMP_ADJUSTED),1,'last').';
for iP = 1 : Data.MetaData.nProfs
    if firstNonNaN(iP,1) > 0
        firstNonNaN(iP,2) = Data.Raw.DEPTH(firstNonNaN(iP,1),iP);
    end
    if lastNonNaN(iP,1) > 0
        lastNonNaN(iP,2) = Data.Raw.DEPTH(lastNonNaN(iP,1),iP);
    end
end

% Exclude profiles from processing where no usable observations are available
ProfileInfo.General.Processed = true(Data.MetaData.nProfs,1);
ProfileInfo.General.Processed(firstNonNaN(:,2) < defaultPars.depthInterpGrid(end)) = false; % no obs above maximum profile depth
ProfileInfo.General.Processed(lastNonNaN(:,2) >= defaultPars.minProfileDepth) = false; % no obs below minimum default profile depth
ProfileInfo.General.Processed(firstNonNaN(:,2) - lastNonNaN(:,2) < defaultPars.minProfileLength) = false; % depth range of usable obs < minimum default profile length

% Distance travelled by seal
distance_km = zeros(Data.MetaData.nProfs,2);
for iP = 2 : Data.MetaData.nProfs
    if all(isfinite([Data.Raw.LATITUDE(iP-1:iP);Data.Raw.LONGITUDE(iP-1:iP)]))
        loc1 = [Data.Raw.LATITUDE(iP-1) Data.Raw.LONGITUDE(iP-1)] ;
        loc2 = [Data.Raw.LATITUDE(iP) Data.Raw.LONGITUDE(iP)] ;
        [dist1km_temp dist2km_temp] = lldistkm(loc1, loc2) ;
        % write computed values in platform_distancekm array
        distance_km(iP,1) = distance_km(iP-1,1) + dist1km_temp ;  % col 1 = cumulated distance
        distance_km(iP,2) = dist1km_temp ;	% col 2 = relative distance
    else
        distance_km(iP,1) = distance_km(iP-1,1) ;
        distance_km(iP,2) = NaN ;
    end
end
ProfileInfo.General.Distance = distance_km(:,2);
ProfileInfo.General.CumDistance = distance_km(:,1);

%% Solar altitude
dateNUM = datenum(dateYMD);
for iP = 1 : Data.MetaData.nProfs
        % compute solar angle and solar altitude (even for all NaN profiles)   
        % solar angle
        zenith_azimuth = solarPosition( ...
            dateNUM(iP),...
            Data.Raw.LATITUDE(iP), ...
            Data.Raw.LONGITUDE(iP),...
            defaultPars.PAR.UTC_offset,0,defaultPars.PAR.DST);

    	% write computed values in genData array
        ProfileInfo.General.SolarAlt(iP) = 90 - zenith_azimuth(1); % solar altitude value (angle, in dregrees)
    	ProfileInfo.General.TimeOfDay_UTC(iP) = datepart(dateYMD(iP),'hour'); % hour of the day (UTC) of the profile
end

% Identify daytime profiles
ProfileInfo.General.Daytime = true(Data.MetaData.nProfs,1);
ProfileInfo.General.Daytime(ProfileInfo.General.SolarAlt <= defaultPars.PAR.SolarAlt_DayTime) = false;

%% Get bathymetry data along seal path
% Extract bathymetry subset around the seal track's lat/lon limits from ETOPO dataset
i_latmin_in_bathy = dsearchn(bathymetry.lat',floor(Data.MetaData.geospatial_lat_min));
i_latmax_in_bathy = dsearchn(bathymetry.lat',ceil(Data.MetaData.geospatial_lat_max));
i_lonmin_in_bathy = dsearchn(bathymetry.lon',floor(convlon(Data.MetaData.geospatial_lon_min,'signed')));
i_lonmax_in_bathy = dsearchn(bathymetry.lon',ceil(convlon(Data.MetaData.geospatial_lon_max,'signed')));

if Data.MetaData.geospatial_lon_max - Data.MetaData.geospatial_lon_min < 180
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
i_lat_in_bathy = dsearchn(bathy_lat',Data.Raw.LATITUDE);
i_lat_in_bathy(~isfinite(Data.Raw.LATITUDE)) = NaN;
i_lon_in_bathy = dsearchn(bathy_lon',Data.Raw.LONGITUDE);
i_lon_in_bathy(~isfinite(Data.Raw.LONGITUDE)) = NaN;
ProfileInfo.General.BathymetryAtProf = NaN(Data.MetaData.nProfs,1);
for iP = 1 : Data.MetaData.nProfs
    if isfinite(i_lat_in_bathy(iP)) && isfinite(i_lon_in_bathy(iP))
        ProfileInfo.General.BathymetryAtProf(iP) = bathy_subset(i_lat_in_bathy(iP),i_lon_in_bathy(iP));
    end
end

%% Get sea ice concentration along seal path
% Load sea ice spatial grid and projection (product: Bremen University ASI-AMSR2 6km)
x = ncread(fullfile(root.data.seaice,'asi-AMSR2-s6250-20210601-v5.4.nc'),'x');
y = ncread(fullfile(root.data.seaice,'asi-AMSR2-s6250-20210601-v5.4.nc'),'y');
proj = projcrs(3412,"Authority","EPSG");
% [X,Y] = meshgrid(x,y);
% [SIC_lat,SIC_lon] = projinv(proj,X,Y);

% Convert profile latitudes and longitudes to x-y coordinates of the sea ice data and find nearest grid point
[x_prof,y_prof] = projfwd(proj,ProfileInfo.General.Lat,ProfileInfo.General.Lon);
i_x_prof = dsearchn(x,x_prof);
i_y_prof = dsearchn(y,y_prof);
% Convert profile date to YYYYMMDD format to find matching sea ice concentration file
profDates = convertTo(ProfileInfo.General.Date,'YYYYMMDD');

% Extract sea ice concentration for each profile at given grid point and day
root_seaice = root.data.seaice;
SeaIceConcAtProf = ProfileInfo.General.SeaIceConcAtProf;
for iP = 1 : numel(profDates)
    SIC_filepath = fullfile(root_seaice,sprintf('asi-AMSR2-s6250-%i-v5.4.nc',profDates(iP)));
    
    if exist(SIC_filepath,'file')
        SIC = ncread(SIC_filepath,'z');
        SeaIceConcAtProf(iP) = SIC(i_x_prof(iP),i_y_prof(iP));
    end
end

ProfileInfo.General.SeaIceConcAtProf = SeaIceConcAtProf;


%% MLD
MLD_algo = 3; % 1: Temperature threshold, 2: Salinity threshold, 3: Density threshold
for iP = 1 : Data.MetaData.nProfs
    % MLD
    pres = Data.Processed.PRES.Reg(:,iP);
    temp = Data.Processed.TEMP.Reg(:,iP);
    psal = Data.Processed.PSAL.Reg(:,iP);
    finiteVals = isfinite(pres) & isfinite(temp) & isfinite(psal);

    % Proceed only if finite T, P, and S values available and if at least five observations
    % are available below the 10 m reference depth (bit of a random number, but normally there should be more data anyway)
    if ~isempty(find(finiteVals,1)) && max(pres(finiteVals)) > 15
        mld_estimates_p = mld(pres(finiteVals),temp(finiteVals),psal(finiteVals),'metric','threshold')';
        mld_estimates = gsw_z_from_p(mld_estimates_p,repmat(ProfileInfo.General.Lat(iP),1,3));
        if mld_estimates(MLD_algo) < 0
            % Use the specified MLD output from the mld function only if depth is <0m
            ProfileInfo.General.MLD(iP) = mld_estimates(MLD_algo);
        else
            % Else use the mean of the T and/or S threshold estimates (whichever yields an MLD <0m)
            ProfileInfo.General.MLD(iP) = mean(mld_estimates(mld_estimates<0));
        end
        
    end
end

%% CMD message: done
pause(0.1)
fprintf('\b\b \x2713\n')

end