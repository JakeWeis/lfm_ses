function [Data,ProfileInfo] = loadData(root,platformID,defaultPars)
% LOADDATA loads raw data and extracts metadata from seal tag or BGC-Argo NetCDF files. Raw data is interpolated onto a
% uniform 1-m depth grid. General profile information is stored in table format.
%
% General profile information:
% 'Profile'                  profile number
% 'PlatformID'               seal tag/float ID
% 'Date'                     profile date
% 'DeployDay'                deployment day
% 'Lon'                      profile longitude
% 'Lat'                      profile latitude
% 'Distance'                 distance travelled by seal/float between profiles
% 'CumDistance'              cumulative distance travelled across all profiles
% 'Processed'                true/false, profile processed (false if no observations available at given profile)
% 'BathymetryAtProf'         bathymetry at profile location
% 'SeaIceConcAtProf'         sea ice concentrations at profile location
% 'TimeOfDay_UTC'            hour of the day (UTC) of the profile
% 'SolarAlt'                 solar altitude value (angle, in dregrees)
% 'Daytime'                  true/false: daytime profile
% 'MLD'                      MLD [T, S, RHO] (based on Holte and Talley 2009 algorithms)
% 'DynamicHeight'            Dynamic height at 50 dbar relative to 1000 dbar (m)
% 'FrontalZone_DynHeight'    Southern Ocean frontal zone based on the dynamic height
% 'FrontalZone_Orsi'         Southern Ocean frontal zone based on Orsi & Harris (2019)
%
% INPUT ARGUMENTS
% root - LFM_SES root directory
%   structure, defined in lfm_ses_menu.m
% platformID - Seal tag or BGC-Argo ID
%   numerical, determined in lfm_ses_menu.m
% defaultPars - default processing parameters
%   structure, defined in setDefaults.m
%
% OUTPUT
% Data - Seal tag/BGC-Argo raw data, processed data and metadata
%   structure, processed data stored after each processing step
% ProfileInfo - Profile specific information/metadata (general, PAR, IRR490, FLUO)
%   structure, profile info listed in table format

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

% Remove bad Argo CTD data (QC flags ~= 1,2,5,8)
qc_keep = [1,2,5,8];
if strcmp(platform_type, 'float')
    Data.Raw.PRES_ADJUSTED(~ismember(str2double(num2cell(Data.Raw.PRES_QC)),qc_keep)) = NaN;
    Data.Raw.TEMP_ADJUSTED(~ismember(str2double(num2cell(Data.Raw.TEMP_QC)),qc_keep)) = NaN;
    Data.Raw.PSAL_ADJUSTED(~ismember(str2double(num2cell(Data.Raw.PSAL_QC)),qc_keep)) = NaN;

    % Remove PAR/F values where 'adjusted' QC flag is 4 (bad)
    Data.Raw.DOWNWELLING_PAR(ismember(str2double(num2cell(Data.Raw.DOWNWELLING_PAR_ADJUSTED_QC)),4)) = NaN;
    Data.Raw.CHLA(ismember(str2double(num2cell(Data.Raw.CHLA_ADJUSTED_QC)),4)) = NaN;
end

% Calculate density from abs. salinity, conservative temperature and pressure (using GSW equations)
Data.Raw.PRES_ADJUSTED(Data.Raw.PRES_ADJUSTED<0) = NaN;   % Set negative pressures (at the surface) NaN
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
        Data.MetaData = ncloadatt_struct(fullfile(root.input, platformID));
        Data.MetaData.platform_type = 'float';
        Data.MetaData.NegativeLight.n_obs_unadj = numel(find(Data.Raw.DOWNWELLING_PAR<0));
        Data.MetaData.NegativeLight.n_prof_unadj = numel(find(any(Data.Raw.DOWNWELLING_PAR<0)));
        Data.MetaData.NegativeLight.n_obs_adj = numel(find(Data.Raw.DOWNWELLING_PAR_ADJUSTED<0));
        Data.MetaData.NegativeLight.n_prof_adj = numel(find(any(Data.Raw.DOWNWELLING_PAR_ADJUSTED<0)));
end

% Number of profiles and number of depth levels
Data.MetaData.nProfs = size(Data.Raw.PRES,2);
Data.MetaData.nDepth = size(Data.Raw.PRES,1);

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
if isfield(Data.Raw,'CHLA')
    Data.Processed.FLUO.Reg = NaN(numel(defaultPars.depthInterpGrid),Data.MetaData.nProfs);
end
if isfield(Data.Raw,'LIGHT')
    Data.Processed.LIGHT.log.Reg = NaN(numel(defaultPars.depthInterpGrid),Data.MetaData.nProfs);
    Data.Processed.LIGHT.lin.Reg = NaN(numel(defaultPars.depthInterpGrid),Data.MetaData.nProfs);
end
if isfield(Data.Raw,'DOWNWELLING_PAR')
    Data.Processed.DOWNWELLING_PAR.log.Reg = NaN(numel(defaultPars.depthInterpGrid),Data.MetaData.nProfs);
    Data.Processed.DOWNWELLING_PAR.lin.Reg = NaN(numel(defaultPars.depthInterpGrid),Data.MetaData.nProfs);
end
if isfield(Data.Raw,'DOWN_IRRADIANCE490')
    Data.Processed.DOWN_IRRADIANCE490.log.Reg = NaN(numel(defaultPars.depthInterpGrid),Data.MetaData.nProfs);
    Data.Processed.DOWN_IRRADIANCE490.lin.Reg = NaN(numel(defaultPars.depthInterpGrid),Data.MetaData.nProfs);
end

% Interpolate
for iP = 1 : Data.MetaData.nProfs
    % depth vector for profile i (irregular step)
    depth_irreg = Data.Raw.DEPTH(:,iP);

    % Identify and remove NaNs and duplicate depth levels
    % --> Duplicate depths can occur in BGC Argo data. 
    iZ_nan = ~isfinite(depth_irreg);
    [~, iUnq] = unique(depth_irreg,'first');
    iZ_dupe = not(ismember(1:numel(depth_irreg),iUnq))';
    depth_irreg(iZ_nan|iZ_dupe) = [];
    
    % Only proceed with interpolation if profile has any finite data
    if numel(depth_irreg)>1
        % PRESSURE
        pres_irreg = Data.Raw.PRES_ADJUSTED(~(iZ_nan|iZ_dupe),iP);
        Data.Processed.PRES.Reg(:,iP) = interp1(depth_irreg,pres_irreg,defaultPars.depthInterpGrid,'linear',NaN);

        % TEMPERATURE
        temp_irreg = Data.Raw.TEMP_ADJUSTED(~(iZ_nan|iZ_dupe),iP);
        Data.Processed.TEMP.Reg(:,iP) = interp1(depth_irreg,temp_irreg,defaultPars.depthInterpGrid,'linear',NaN);

        % PRACTICAL SALINITY
        psal_irreg = Data.Raw.PSAL_ADJUSTED(~(iZ_nan|iZ_dupe),iP);
        Data.Processed.PSAL.Reg(:,iP) = interp1(depth_irreg,psal_irreg,defaultPars.depthInterpGrid,'linear',NaN);

        % CHL FLUORESCENCE
        if isfield(Data.Raw,'CHLA')
            fluo_irreg = Data.Raw.CHLA(~(iZ_nan|iZ_dupe),iP);
            if strcmp(platform_type, 'sealtag')
                Data.Processed.FLUO.Reg(:,iP) = interp1(depth_irreg,fluo_irreg,defaultPars.depthInterpGrid,'linear',NaN);
            elseif strcmp(platform_type, 'float')
                % Only consider depths at which Fluo observations are available 
                % (different to seal tag data, BGC variables are not measured at every available depth level)
                fluo_isnan = ~isfinite(fluo_irreg);
                if numel(find(~fluo_isnan))>1
                    Data.Processed.FLUO.Reg(:,iP) = interp1(depth_irreg(~fluo_isnan),movmedian(fluo_irreg(~fluo_isnan),5,'omitnan'),defaultPars.depthInterpGrid,'linear',NaN);
                end
            end
        end

        % LIGHT/PAR
        % (1) Seal tag (light)
        % NB: log-transformed, linear interpolation = exponential decrease to depth
        if strcmp(platform_type, 'sealtag') && isfield(Data.Raw,'LIGHT')
            % Seal tag light data is ln(light)
            parlog_irreg = Data.Raw.LIGHT(~(iZ_nan|iZ_dupe),iP);
            Data.Processed.LIGHT.log.Reg(:,iP) = interp1(depth_irreg,parlog_irreg,defaultPars.depthInterpGrid,'linear',NaN);
            % Reverse log-transformation
            Data.Processed.LIGHT.lin.Reg(:,iP) = exp(Data.Processed.LIGHT.log.Reg(:,iP));
        end

        % (2): BGC Argo
        % NB: needs to be log-transformed prior to linear interpolation
        if strcmp(platform_type, 'float')
            % Downwelling PAR
            if isfield(Data.Raw,'DOWNWELLING_PAR')
                par_irreg = Data.Raw.DOWNWELLING_PAR(~(iZ_nan|iZ_dupe),iP);
                par_irreg(par_irreg<=0) = NaN;      % remove negative/0 values before log-transformation
                parlog_irreg = log(par_irreg);
                
                % Only consider depths at which PAR observations are available
                % NB: Float BGC variables are not measured at every available depth level 
                i_par_isnan = ~isfinite(parlog_irreg);
                if numel(find(~i_par_isnan))>1
                    Data.Processed.DOWNWELLING_PAR.log.Reg(:,iP) = interp1(depth_irreg(~i_par_isnan),parlog_irreg(~i_par_isnan),defaultPars.depthInterpGrid,'linear',NaN);
                    Data.Processed.DOWNWELLING_PAR.lin.Reg(:,iP) = exp(Data.Processed.DOWNWELLING_PAR.log.Reg(:,iP));
                end
            end
            
            % BGC-Argo downwelling irradiance at 490 nm
            if isfield(Data.Raw,'DOWN_IRRADIANCE490')
                irr490_irreg = Data.Raw.DOWN_IRRADIANCE490(~(iZ_nan|iZ_dupe),iP);
                irr490_irreg(irr490_irreg<=0) = NaN;      % remove negative/0 values before log-transformation
                irr490log_irreg = log(irr490_irreg);
                
                % Only consider depths at which PAR observations are available
                % NB: Float BGC variables are not measured at every available depth level 
                i_par_isnan = ~isfinite(irr490log_irreg);
                if numel(find(~i_par_isnan))>1
                    Data.Processed.DOWN_IRRADIANCE490.log.Reg(:,iP) = interp1(depth_irreg(~i_par_isnan),irr490log_irreg(~i_par_isnan),defaultPars.depthInterpGrid,'linear',NaN);
                    Data.Processed.DOWN_IRRADIANCE490.lin.Reg(:,iP) = exp(Data.Processed.DOWN_IRRADIANCE490.log.Reg(:,iP));
                end
            end
        end
    end
end

%% create genData table to store general profile information
var_names = {...
    'Profile', ...                  % profile number
    'PlatformID', ...               % seal tag/float ID
    'Date', ...                     % date
    'DeployDay', ...                % deployment day
    'Lon', ...                      % lon
    'Lat', ...                      % lat
    'Distance', ...                 % distance travelled by seal between profiles
    'CumDistance', ...              % cumulative distance travelled by seal across all profiles
    'Processed', ...                % true/false: profile processed
    'Bathymetry_BEDMAP', ...        % bathymetry from BEDMAP2
    'Bathymetry_ETOPO', ...         % bathymetry from ETOPO 2022 (0.5˚/30s res.)
    'SeaIceConc', ...               % sea ice concentrations at profile location
    'SeaIceCovered', ...            % true/false if profile was potentialy ice covered
    'TimeOfDay_UTC', ...            % hour of the day (UTC) of the profile
    'SolarAlt', ...                 % solar altitude value (angle, in dregrees)
    'Daytime', ...                  % true/false: daytime profile
    'MLD', ...                      % MLD [T, S, RHO] (based on Holte and Talley 2009 algorithms)
    'DynamicHeight', ...            % Dynamic height at 50 dbar relative to 1000 dbar (m)
    'FrontalZone_DynHeight',...     % Southern Ocean frontal zone based on the dynamic height
    'FrontalZone_Orsi'              % Southern Ocean frontal zone based on Orsi & Harris (2019)
    }';

ProfileInfo.General = array2table(NaN(Data.MetaData.nProfs, numel(var_names)),'VariableNames',var_names);

% Write data to genData table
ProfileInfo.General.Profile = (1 : Data.MetaData.nProfs)';
switch platform_type
    case 'sealtag'
        ProfileInfo.General.PlatformID = repmat(Data.MetaData.smru_platform_code,Data.MetaData.nProfs,1);
    case 'float'
        ProfileInfo.General.PlatformID = ncread(fullfile(root.input, platformID),'PLATFORM_NUMBER')';
end
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
dateNUM = convertTo(dateYMD,'datenum');
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

%% Get bathymetry data along seal path from ETOPO 2022 and BEDMAP2
missingBathyFiles = (...
    ~exist(fullfile(root.data.etopo,'ETOPO_2022_v1_30s_N90W180_bed.nc'),"file") &&...
    ~exist(fullfile(root.data.etopo,'ETOPO_2022_v1_0s_N90W180_bed.nc'),"file")) ||...
    ~exist(fullfile(root.data.bedmap,'bedmap2_bed.flt'),"file") ||...
    ~exist(fullfile(root.data.bedmap,'bedmap2_grounded_bed_uncertainty.flt'),"file");

if ~missingBathyFiles
    [ProfileInfo.General.Bathymetry_ETOPO,ProfileInfo.General.Bathymetry_BEDMAP] = ...
        ll2bathy(ProfileInfo.General.Lat,ProfileInfo.General.Lon,root);
end

%% Get sea ice concentration along seal path
missingSICFiles = isempty(dirPaths(fullfile(root.data.seaice,'*.nc')));
if ~missingSICFiles
    allFiles = dirPaths([root.data.seaice,filesep,'*.nc']);
    % Load sea ice spatial grid and projection (product: Bremen University ASI-AMSR2 6km)
    x = ncread(allFiles(1).path,'x');
    y = ncread(allFiles(1).path,'y');
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
    SeaIceConcAtProf = ProfileInfo.General.SeaIceConc;
    for iP = 1 : numel(profDates)
        SIC_filepath = fullfile(root_seaice,sprintf('asi-AMSR2-s6250-%i-v5.4.nc',profDates(iP)));

        if exist(SIC_filepath,'file')
            SIC = ncread(SIC_filepath,'z');
            SeaIceConcAtProf(iP) = SIC(i_x_prof(iP),i_y_prof(iP));
        end
    end

    ProfileInfo.General.SeaIceConc = SeaIceConcAtProf;
end

%% Sea ice detection
% Following Riser, S. C., Swift, D., and Drucker, R., 2018 (https://doi.org/10.1002/2017JC013419)
% The ice-sensing algorithm simply compares the median temperature between ∼ 50 and 20 m during ascent to a threshold temperature of −1.78 ∘C.
ProfileInfo.General.SeaIceCovered = median(Data.Processed.TEMP.Reg(2:50,:))' <= defaultPars.SI_T_threshold;

%% MLD, dynamic height & SO frontal zone
MLD_algo = 3; % 1: Temperature threshold, 2: Salinity threshold, 3: Density threshold
for iP = 1 : Data.MetaData.nProfs
    % MLD
    pres = Data.Processed.PRES.Reg(:,iP);
    temp = Data.Processed.TEMP.Reg(:,iP);
    psal = Data.Processed.PSAL.Reg(:,iP);
    finiteVals = isfinite(pres) & isfinite(temp) & isfinite(psal);

    % Proceed only if finite T, P, and S values available and if at least 3 observations are available below the 10 m reference depth
    if numel(find(pres(finiteVals)>10))>3
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

    % Dynamic height
    lon = ProfileInfo.General.Lon(iP);
    lat = ProfileInfo.General.Lat(iP);
    asal = gsw_SA_from_SP(psal, pres, lon, lat);
    ct = gsw_CT_from_t(asal,temp,pres);

    finiteVals = isfinite(pres) & isfinite(temp) & isfinite(psal);
    pres_finiteObs = pres;
    pres_finiteObs(~finiteVals) = NaN;

    iZ_50dbar = dsearchn(pres_finiteObs,50);
    iZ_1000dbar = dsearchn(pres_finiteObs,1000);
    
    % Only proceed if the target (50 dbar) and reference pressure (1000 dbar) are close enough to literature values
    if iZ_50dbar < 100 && iZ_1000dbar > 950
        g = 9.81;   % Gravitational constant
        dynHeight = gsw_geo_strf_dyn_height(asal(iZ_50dbar:iZ_1000dbar),ct(iZ_50dbar:iZ_1000dbar),pres(iZ_50dbar:iZ_1000dbar),pres(iZ_1000dbar))/g;
        ProfileInfo.General.DynamicHeight(iP) = dynHeight(1);

        % SO frontal zones
        % SAF   --> 0.87 dynamic height
        % PF    --> 0.62 dynamic height
        % sACCf --> 0.45 dynamic height
        % sBdy  --> 0.38 dynamic height
        if ProfileInfo.General.DynamicHeight(iP) >= 0.87
            ProfileInfo.General.FrontalZone_DynHeight(iP) = 1;
        elseif ProfileInfo.General.DynamicHeight(iP) < 0.87 && ProfileInfo.General.DynamicHeight(iP) >= 0.62
            ProfileInfo.General.FrontalZone_DynHeight(iP) = 2;
        elseif ProfileInfo.General.DynamicHeight(iP) < 0.62 && ProfileInfo.General.DynamicHeight(iP) >= 0.45
            ProfileInfo.General.FrontalZone_DynHeight(iP) = 3;
        elseif ProfileInfo.General.DynamicHeight(iP) < 0.45 && ProfileInfo.General.DynamicHeight(iP) >= 0.38
            ProfileInfo.General.FrontalZone_DynHeight(iP) = 4;
        elseif ProfileInfo.General.DynamicHeight(iP) < 0.38
            ProfileInfo.General.FrontalZone_DynHeight(iP) = 5;
        else
            ProfileInfo.General.FrontalZone_DynHeight(iP) = NaN;
        end
    end
end

% Determined fronts based on average front and zone locations (Orsi & Harris, 2019)
[~,zones] = SOFronts;
ProfileInfo.General.FrontalZone_Orsi(inpolygon(ProfileInfo.General.Lon,ProfileInfo.General.Lat,zones.SAZ.lon,zones.SAZ.lat)) = 1;
ProfileInfo.General.FrontalZone_Orsi(inpolygon(ProfileInfo.General.Lon,ProfileInfo.General.Lat,zones.PFZ.lon,zones.PFZ.lat)) = 2;
ProfileInfo.General.FrontalZone_Orsi(inpolygon(ProfileInfo.General.Lon,ProfileInfo.General.Lat,zones.AZ.lon,zones.AZ.lat)) = 3;
ProfileInfo.General.FrontalZone_Orsi(inpolygon(ProfileInfo.General.Lon,ProfileInfo.General.Lat,zones.SZ.lon,zones.SZ.lat)) = 4;
ProfileInfo.General.FrontalZone_Orsi(inpolygon(ProfileInfo.General.Lon,ProfileInfo.General.Lat,zones.SPR.lon,zones.SPR.lat)) = 5;

ProfileInfo.General.FrontalZone_DynHeight = categorical(ProfileInfo.General.FrontalZone_DynHeight,1:5,{'SAZ','PFZ','AZ','SZ','SPR'});
ProfileInfo.General.FrontalZone_Orsi = categorical(ProfileInfo.General.FrontalZone_Orsi,1:5,{'SAZ','PFZ','AZ','SZ','SPR'});


%% CMD message: done
pause(0.1)
fprintf('\b\b \x2713\n')


%% Custom functions
function [B_etopo,B_bedmap,B_unc_bedmap] = ll2bathy(lat,lon,root)
% LLBATHY extracts the bathymetry depth at given coordinates.
%
% Input: lat and lon coordinates
% 
% Optional: 
% m_ll2bathy(lat, lon, XDATA), plot bathymetric data against given XDATA (e.g., date, lat, lon). Must be the same length of
% lat/lon data).
%
% Output:
% B_bedmap - bathymetry based on BEDMAP2 (NaN outside polarstereographic grid limits)
% B_etopo - bathymetry based on ETOPO 2022

% Reshape lat lon to be column vectors
assert(isvector(lat)&isvector(lon),'Lat and lon must be vectors.')
lat = reshape(lat,numel(lat),1);
lon = reshape(lon,numel(lon),1);

%% BEDMAP2 BATHYMETRY
% Convert query lat/lon coordinates to polar steroegraphic x/y coordinates
x_query = ones(numel(lat),1)*3333500+1;
y_query = ones(numel(lat),1)*3333500+1;
for iL = 1 : numel(lat)
    if lat(iL)<0
        % ll2ps throws an error if the latitude is in the northern hemisphere.
        [x_query(iL),y_query(iL)] = ll2ps(lat(iL),lon(iL));
    end
end

% Polarstereographic grid
x_bedmap = linspace(-3333500,3333500,6667)';
y_bedmap = linspace(-3333500,3333500,6667)';

% Open binary bathymetry and uncertainty files (.flt)
fid_1 = fopen(fullfile(root.data.bedmap,'bedmap2_bed.flt'), 'r', 'l');
fid_2 = fopen(fullfile(root.data.bedmap,'bedmap2_grounded_bed_uncertainty.flt'), 'r', 'l');

% row/column indices of the query points (x/y PS-grid points calculated from lat/lon coordinates) 
% NB: The bathymetry grid in the .flt file is rotated by 90 degrees (CW) relative to the polar stereographic grid. Therefore, finding the
% correct point in the file requires going "backwards" along the columns. This is why the column index is subtracted from the
% total number of columns (6667). The row index does not change with the rotation.
i_row = dsearchn(x_bedmap,x_query);
i_col = 6667-dsearchn(y_bedmap,y_query)+1;

% Offset required to get to the query point in the file (moving along columns).
offset = (i_col - 1) * 6667 + (i_row - 1);

% Read bathymetry and uncertainty values from files
B_bedmap = NaN(numel(i_row),1);
B_unc_bedmap = NaN(numel(i_row),1);
for iL = 1 : numel(i_row)
    if (x_query(iL)>=-3333500 && x_query(iL)<=3333500) && (y_query(iL)>=-3333500 && y_query(iL)<=3333500)
        fseek(fid_1,offset(iL)*4,'bof');
        B_bedmap(iL) = fread(fid_1, [1 1], 'float32');
        fseek(fid_2,offset(iL)*4,'bof');
        B_unc_bedmap(iL) = fread(fid_2, [1 1], 'float32');
    else
        B_bedmap(iL) = NaN;
        B_unc_bedmap(iL) = NaN;
    end
end
fclose(fid_1);
fclose(fid_2);

%% ETOPO 2022 BATHYMETRY

% Extract data from netCDF file at specified coordinates
lat_etopo = ncread(fullfile(root.data.etopo,'ETOPO_2022_v1_30s_N90W180_bed.nc'),'lat');
lon_etopo = ncread(fullfile(root.data.etopo,'ETOPO_2022_v1_30s_N90W180_bed.nc'),'lon');
i_lat = dsearchn(lat_etopo,lat);
i_lon = dsearchn(lon_etopo,lon);

B_etopo = NaN(numel(i_lat),1);
for iL = 1 : numel(i_lat)
    B_etopo(iL) = ncread(fullfile(root.data.etopo,'ETOPO_2022_v1_30s_N90W180_bed.nc'),'z',[i_lon(iL),i_lat(iL)],[1 1]);
end

end

function [Fronts,Zones] = SOFronts(varargin)
% Extract Orsi, AH., Harris, U. (2019) Southern Ocean front coordinates. Coordinates are centered around the 0˚ meridian.
%
% [fronts,zones] = SOFronts;
% Returns standard Orsi fronts and zones:
% >> NBdy       Northern Southern Ocean boundary (30˚S)
%    STZ        Subtropical zone
% >> STF        Subtropical front
%    SAZ        Subantarctic Zone
% >> SAF        Subantarctic front
%    PFZ        Polar frontal zone
% >> PF         Polar front
%    AZ         Antarctic Zone
% >> SACCF      Southern Antarctic Circumpolar Current front
%    SZ         Southern zone
% >> SBdy       Southern boundary
%    SPR        Subpolar region
% 
% Add custom fronts (multiple fronts as n-by-3 cell array):
% 'add_fronts'          {'front_name',longitudes,latitudes}
%   front_name          Name of the custom front (NB: used as structure field name)
%   longitudes          Longitudes of the custom front
%   latitudes           Latitudes of the custom front
% 
% Add custom zone (multiple zones as n-by-3 cell array)
% 'add_zone'            {'zone_name','N_front','S_front'}
%   zone_name           Name of the custom zone  (NB: used as structure field name)
%   N_front             Name of the northern front  (either one of the standard Orsi fronts or a custom front)
%   S_front             Name of the southern front  (either one of the standard Orsi fronts or a custom front)
%
% Example:
% front1 = {'lat_53S',linspace(-180,180,500),repmat(-53,500,1)};
% front2 = {'lat_63S',linspace(-180,180,500),repmat(-63,500,1)};
% custom_fronts = [front1;front2];
% zone1 = {'STF_53S','STF','lat_53S'};
% zone2 = {'SAF_63S','SAF','lat_63S'};
% custom_zones = [zone1,zone2];
%
% [fronts,zones] = SOFronts('add_fronts',custom_fronts,'add_zones',custom_zones);
%
% Reference: Orsi, AH., Harris, U. (2019) Fronts of the Antarctic Circumpolar Current GIS data, Ver. 1, Australian Antarctic
% Data Centre https://data.aad.gov.au/metadata/records/antarctic_circumpolar_current_fronts Accessed: 2020-11-19

% Parse input
p = inputParser;
addParameter(p,'add_fronts',{})
addParameter(p,'add_zones',{})
parse(p,varargin{:});

opt.AddFronts = p.Results.add_fronts;
opt.AddZones = p.Results.add_zones;

frontstable = readtable(fullfile(root.proj,'/02_library/OrsiFronts.csv'));

% SO front lat/lon data
front_name = {'NBdy','STF','SAF','PF','SACCF','SBdy'};
for iF = 1 : length(front_name)
    % Extract from Orsi et al data table
    Fronts.(front_name{iF}).lon = frontstable{strcmp(frontstable{:,3},front_name{iF}),1};
    Fronts.(front_name{iF}).lat = frontstable{strcmp(frontstable{:,3},front_name{iF}),2};

    % Extend vectors to -180 and 180˚ lon (repeating first and last latitude value)
    Fronts.(front_name{iF}).lon = [-180;Fronts.(front_name{iF}).lon;180];
    Fronts.(front_name{iF}).lat = [Fronts.(front_name{iF}).lat(1);Fronts.(front_name{iF}).lat;Fronts.(front_name{iF}).lat(end)];
end

% SO frontal zone polygons
zone_name = {'STZ','SAZ','PFZ','AZ','SZ','SPR'};
for iZ = 1 : numel(zone_name)
    if iZ < numel(zone_name)
        Zones.(zone_name{iZ}).lon = [Fronts.(front_name{iZ}).lon;flipud(Fronts.(front_name{iZ+1}).lon)];
        Zones.(zone_name{iZ}).lat = [Fronts.(front_name{iZ}).lat;flipud(Fronts.(front_name{iZ+1}).lat)];
    else
        Zones.(zone_name{iZ}).lon = [Fronts.(front_name{iZ}).lon;linspace(180,-180,100)'];
        Zones.(zone_name{iZ}).lat = [Fronts.(front_name{iZ}).lat;repmat(-89,100,1)];
    end
end

% Add custom fronts
if ~isempty(opt.AddFronts)
    for iCF = 1 : size(opt.AddFronts,1)
        Fronts.(opt.AddFronts{iCF,1}).lon = reshape(opt.AddFronts{iCF,2},numel(opt.AddFronts{iCF,2}),1);
        Fronts.(opt.AddFronts{iCF,1}).lat = reshape(opt.AddFronts{iCF,3},numel(opt.AddFronts{iCF,3}),1);
    end
end

% Add custom zones
if ~isempty(opt.AddZones)
    for iCZ = 1 : size(opt.AddZones,1)
        Zones.(opt.AddZones{iCZ,1}).lon = [Fronts.(opt.AddZones{iCZ,2}).lon;flipud(Fronts.(opt.AddZones{iCZ,3}).lon)];
        Zones.(opt.AddZones{iCZ,1}).lat = [Fronts.(opt.AddZones{iCZ,2}).lat;flipud(Fronts.(opt.AddZones{iCZ,3}).lat)];
    end
end

end

end