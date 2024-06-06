function [tagData,tagMetadata,tagProcessed,ProfileInfo] = loadData(root,tagRef,bathymetry,defaultPars)
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
if strcmp(ncreadatt(fullfile(root.input, tagRef),'/','source'), 'Marine mammal observation')
    platform_type = 'sealtag';
elseif strcmp(ncreadatt(fullfile(root.input, tagRef),'/','source'), 'Argo float')
    platform_type = 'float';
else
    platform_type = 'undefined';
    warning('File is neither seal tag nor Argo float.')
end

%% Load data
tagData = ncload_struct(fullfile(root.input, tagRef));

% If available, replace LATITUDE/LONGITUDE fields with corresponding *_ADJUSTED fields
if isfield(tagData,'LATITUDE_ADJUSTED')
    tagData.LATITUDE = tagData.LATITUDE_ADJUSTED;
    tagData.LONGITUDE = tagData.LONGITUDE_ADJUSTED;
end

% Convert JULD to DATETIME and DATENUM
reference_date = '1950-01-01 00:00:00'; % see platform_data.REFERENCE_DATE_TIME
dateYMD = datetime(tagData.JULD_LOCATION + datenum(reference_date),'ConvertFrom','datenum');
tagData.DATETIME = repmat(dateYMD',size(tagData.PRES,1),1);

% Calculate density from abs. salinity, conservative temperature and pressure (using GSW equations)
SA = gsw_SA_from_SP(tagData.PSAL_ADJUSTED,tagData.PRES_ADJUSTED,tagData.LONGITUDE,tagData.LATITUDE);
CT = gsw_CT_from_t(SA,tagData.TEMP_ADJUSTED,tagData.PRES_ADJUSTED);
tagData.DENS = gsw_rho_CT_exact(SA,CT,0);

%% Load metadata
switch platform_type
    case 'sealtag'
        % Seal tag metadata stored in main data file
        tagMetadata = ncloadatt_struct(fullfile(root.input, tagRef));
        tagMetadata.platform_type = 'sealtag';
    case 'float'
        % Float metadata stored in separate file
        float_metadata_filename = strrep(tagRef,'Sprof','meta');
        tagMetadata = ncload_struct(fullfile(root.input, tagRef));
        tagMetadata.platform_type = 'float';
end

% Number of profiles and number of depth levels
tagMetadata.nProfs = size(tagData.PRES,2);
tagMetadata.nDepth = size(tagData.PRES,1);

% Check resolution of tag data
if ~isempty(regexp(tagRef,'lr[1-9]','once'))
    % Low res: If the file name contains "lr1/lr2/..."
    tagMetadata.depth_res = 'lr';
elseif ~isempty(regexp(tagRef,'hr[1-9]','once'))
    % High res: If the file name contains "hr1/hr2/..."
    tagMetadata.depth_res = 'hr';
elseif ~isempty(regexp(tagRef,'fr[1-9]','once'))
    % Full Res: If the file name contains "fr1/fr2/..."
    tagMetadata.depth_res = 'fr';
else
    % If the resolution is not specified in the file name.
    if tagMetadata.nDepth < numel(defaultPars.depthInterpGrid)
        % Low res: if n depth levels < n depth levels of the specified depth interpolation vector
        tagMetadata.depth_res = 'lr';
    elseif tagMetadata.nDepth == numel(defaultPars.depthInterpGrid)
        % if n depth levels = n depth levels of the specified depth interpolation vector
        tagMetadata.depth_res = 'hr';
    end
end

% lat/lon min/max
tagMetadata.geospatial_lat_max = max(tagData.LATITUDE);
tagMetadata.geospatial_lat_min = min(tagData.LATITUDE);
tagMetadata.geospatial_lon_max = max(tagData.LONGITUDE);
tagMetadata.geospatial_lon_min = min(tagData.LONGITUDE);

% sort depth from deepest to shallowest and find the deepest depth where temperature is not NaN
PRES_vec = tagData.PRES(:);
[~,i_sort] = sort(PRES_vec,'descend','MissingPlacement','last');
TEMP_vec = tagData.TEMP(:);
tagMetadata.max_depth_dataset = max(PRES_vec);
for iP = 1 : numel(i_sort)
    tagMetadata.max_depth_dataset = PRES_vec(i_sort(iP));
    if isfinite(TEMP_vec(i_sort(iP)))
        break
    end
end

%% Interpolate data onto uniform depth grid

% Calculate depth from pressure and latitude (using GSW equations)
tagData.DEPTH = gsw_z_from_p(tagData.PRES_ADJUSTED,tagData.LATITUDE);

% The depth array of the processed data is equivalent to the predefined interpolation-depth vector (in meters)
tagProcessed.DEPTH = repmat(defaultPars.depthInterpGrid,1,tagMetadata.nProfs);

% Pressure, temperature, salinity, fluorescence and PAR are all interpolated onto the uniform depth grid given above
tagProcessed.PRES.Reg = NaN(numel(defaultPars.depthInterpGrid),tagMetadata.nProfs);
tagProcessed.TEMP.Reg = NaN(numel(defaultPars.depthInterpGrid),tagMetadata.nProfs);
tagProcessed.PSAL.Reg = NaN(numel(defaultPars.depthInterpGrid),tagMetadata.nProfs);
tagProcessed.FLUO.Reg = NaN(numel(defaultPars.depthInterpGrid),tagMetadata.nProfs);
tagProcessed.PAR.log.Reg = NaN(numel(defaultPars.depthInterpGrid),tagMetadata.nProfs);
tagProcessed.PAR.lin.Reg = NaN(numel(defaultPars.depthInterpGrid),tagMetadata.nProfs);

% Interpolate
for iP = 1 : tagMetadata.nProfs
    % depth vector for profile i (irregular step)
    depth_irreg = tagData.DEPTH(:,iP);
    depth_isnan = ~isfinite(depth_irreg);
    depth_irreg(depth_isnan) = [];

    % Pressure
    pres_irreg = tagData.PRES_ADJUSTED(:,iP);
    pres_irreg(depth_isnan) = [];
    tagProcessed.PRES.Reg(:,iP) = interp1(depth_irreg,pres_irreg,defaultPars.depthInterpGrid,'linear');

    % Temperature
    temp_irreg = tagData.TEMP_ADJUSTED(:,iP);
    temp_irreg(depth_isnan) = [];
    tagProcessed.TEMP.Reg(:,iP) = interp1(depth_irreg,temp_irreg,defaultPars.depthInterpGrid,'linear');

    % Practical salinity
    psal_irreg = tagData.PSAL_ADJUSTED(:,iP);
    psal_irreg(depth_isnan) = [];
    tagProcessed.PSAL.Reg(:,iP) = interp1(depth_irreg,psal_irreg,defaultPars.depthInterpGrid,'linear');

    % Chl fluorescence
    if isfield(tagData,'CHLA')
        fluo_irreg = tagData.CHLA(:,iP);
        fluo_irreg(depth_isnan) = [];
        tagProcessed.FLUO.Reg(:,iP) = interp1(depth_irreg,fluo_irreg,defaultPars.depthInterpGrid,'linear');
    end

    % LIGHT fluorescence
    if isfield(tagData,'LIGHT')
        par_irreg = tagData.LIGHT(:,iP);
        par_irreg(depth_isnan) = [];
        tagProcessed.PAR.log.Reg(:,iP) = interp1(depth_irreg,par_irreg,defaultPars.depthInterpGrid,'linear');
        tagProcessed.PAR.lin.Reg(:,iP) = exp(tagProcessed.PAR.log.Reg(:,iP));
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
    'TimeOfDay_UTC', ...            % hour of the day (UTC) of the profile
    'SolarAlt', ...                 % solar altitude value (angle, in dregrees)
    'Daytime', ...                  % true/false: daytime profile
    'MLD', ...                      % MLD [T, S, RHO] (based on Holte and Talley 2009 algorithms)
    }';

ProfileInfo = array2table(NaN(tagMetadata.nProfs, numel(var_names)),'VariableNames',var_names);

% Write data to genData table
ProfileInfo.Profile = (1 : tagMetadata.nProfs)';
ProfileInfo.TagID = repmat(tagMetadata.smru_platform_code,tagMetadata.nProfs,1);
ProfileInfo.Date = dateYMD;
ProfileInfo.DeployDay = floor(datenum(dateYMD)) - floor(datenum(dateYMD(1))) + 1;
ProfileInfo.Lon = tagData.LONGITUDE;
ProfileInfo.Lat = tagData.LATITUDE;

% Shallowest and deepest available observation (based on temperature)
firstNonNaN = NaN(tagMetadata.nProfs,2);
lastNonNaN = NaN(tagMetadata.nProfs,2);
firstNonNaN(:,1) = find_ndim(isfinite(tagData.TEMP_ADJUSTED),1,'first').';
lastNonNaN(:,1) = find_ndim(isfinite(tagData.TEMP_ADJUSTED),1,'last').';
for iP = 1 : tagMetadata.nProfs
    if firstNonNaN(iP,1) > 0
        firstNonNaN(iP,2) = tagData.DEPTH(firstNonNaN(iP,1),iP);
    end
    if lastNonNaN(iP,1) > 0
        lastNonNaN(iP,2) = tagData.DEPTH(lastNonNaN(iP,1),iP);
    end
end

% Exclude profiles from processing where no usable observations are available
ProfileInfo.Processed = true(tagMetadata.nProfs,1);
ProfileInfo.Processed(firstNonNaN(:,2) < defaultPars.depthInterpGrid(end)) = false; % no obs above maximum profile depth
ProfileInfo.Processed(lastNonNaN(:,2) >= defaultPars.minProfileDepth) = false; % no obs below minimum default profile depth
ProfileInfo.Processed(firstNonNaN(:,2) - lastNonNaN(:,2) < defaultPars.minProfileLength) = false; % depth range of usable obs < minimum default profile length

% Distance travelled by seal
distance_km = zeros(tagMetadata.nProfs,2);
for iP = 2 : tagMetadata.nProfs
    if all(isfinite([tagData.LATITUDE(iP-1:iP);tagData.LONGITUDE(iP-1:iP)]))
        loc1 = [tagData.LATITUDE(iP-1) tagData.LONGITUDE(iP-1)] ;
        loc2 = [tagData.LATITUDE(iP) tagData.LONGITUDE(iP)] ;
        [dist1km_temp dist2km_temp] = lldistkm(loc1, loc2) ;
        % write computed values in platform_distancekm array
        distance_km(iP,1) = distance_km(iP-1,1) + dist1km_temp ;  % col 1 = cumulated distance
        distance_km(iP,2) = dist1km_temp ;	% col 2 = relative distance
    else
        distance_km(iP,1) = distance_km(iP-1,1) ;
        distance_km(iP,2) = NaN ;
    end
end
ProfileInfo.Distance = distance_km(:,2);
ProfileInfo.CumDistance = distance_km(:,1);

%% Solar altitude
dateNUM = datenum(dateYMD);
for iP = 1 : tagMetadata.nProfs
        % compute solar angle and solar altitude (even for all NaN profiles)   
        % solar angle
        zenith_azimuth = solarPosition( ...
            dateNUM(iP),...
            tagData.LATITUDE(iP), ...
            tagData.LONGITUDE(iP),...
            defaultPars.PAR.UTC_offset,0,defaultPars.PAR.DST);

    	% write computed values in genData array
        ProfileInfo.SolarAlt(iP) = 90 - zenith_azimuth(1); % solar altitude value (angle, in dregrees)
    	ProfileInfo.TimeOfDay_UTC(iP) = datepart(dateYMD(iP),'hour'); % hour of the day (UTC) of the profile
end

% Identify daytime profiles
ProfileInfo.Daytime = true(tagMetadata.nProfs,1);
ProfileInfo.Daytime(ProfileInfo.SolarAlt <= defaultPars.PAR.SolarAlt_DayTime) = false;

%% Get bathymetry data along seal path
% Extract bathymetry subset around the seal track's lat/lon limits from ETOPO dataset
i_latmin_in_bathy = dsearchn(bathymetry.lat',floor(tagMetadata.geospatial_lat_min));
i_latmax_in_bathy = dsearchn(bathymetry.lat',ceil(tagMetadata.geospatial_lat_max));
i_lonmin_in_bathy = dsearchn(bathymetry.lon',floor(convlon(tagMetadata.geospatial_lon_min,'signed')));
i_lonmax_in_bathy = dsearchn(bathymetry.lon',ceil(convlon(tagMetadata.geospatial_lon_max,'signed')));

if tagMetadata.geospatial_lon_max - tagMetadata.geospatial_lon_min < 180
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
i_lat_in_bathy = dsearchn(bathy_lat',tagData.LATITUDE);
i_lat_in_bathy(~isfinite(tagData.LATITUDE)) = NaN;
i_lon_in_bathy = dsearchn(bathy_lon',tagData.LONGITUDE);
i_lon_in_bathy(~isfinite(tagData.LONGITUDE)) = NaN;
ProfileInfo.BathymetryAtProf = NaN(tagMetadata.nProfs,1);
for iP = 1 : tagMetadata.nProfs
    if isfinite(i_lat_in_bathy(iP)) && isfinite(i_lon_in_bathy(iP))
        ProfileInfo.BathymetryAtProf(iP) = bathy_subset(i_lat_in_bathy(iP),i_lon_in_bathy(iP));
    end
end

%% MLD
MLD_algo = 3; % 1: Temperature threshold, 2: Salinity threshold, 3: Density threshold
for iP = 1 : tagMetadata.nProfs
    % MLD
    pres = tagProcessed.PRES.Reg(:,iP);
    temp = tagProcessed.TEMP.Reg(:,iP);
    psal = tagProcessed.PSAL.Reg(:,iP);
    finiteVals = isfinite(pres) & isfinite(temp) & isfinite(psal);

    % Proceed only if finite T, P, and S values available and if at least one observation 
    % is available beloew the 10 m reference depth 
    if ~isempty(find(finiteVals,1)) && max(pres(finiteVals)) > 10
        mld_estimates_p = mld(pres(finiteVals),temp(finiteVals),psal(finiteVals),'metric','threshold')';
        mld_estimates = gsw_z_from_p(mld_estimates_p,repmat(ProfileInfo.Lat(iP),1,3));
        if mld_estimates(MLD_algo) < 0
            % Use the specified MLD output from the mld function only if depth is <0m
            ProfileInfo.MLD(iP) = mld_estimates(MLD_algo);
        else
            % Else use the mean of the T and/or S threshold estimates (whichever yields an MLD <0m)
            ProfileInfo.MLD(iP) = mean(mld_estimates(mld_estimates<0));
        end
        
    end
end

%% CMD message: done
pause(0.1)
fprintf('\b\b \x2713\n')

end