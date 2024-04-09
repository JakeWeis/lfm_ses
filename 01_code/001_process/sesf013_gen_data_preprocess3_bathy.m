%% script metadata
%
% SCRIPT TEMPORARILY DISABLED (SEE CHANGE LOG FOR FURTHER INFORMATION)
%
% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% compute bathymetry along platform trajectory
% when LAT/LON is available
% highlight profiles where maximum depth recorded by platform > bathymetry data
% store in platform_metadata structure

%%%%%%%%%% REFERENCES %%%%%%%%%%
% ETOPO1
% https://www.ngdc.noaa.gov/mgg/global/
% Cite ETOPO1: doi:10.7289/V5C8276M

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x
% sesf010_set_date_pressure_arrays

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])

% +++ CHANGE LOG +++
% L. Le Ster (lls)
% J. Weis (jw)
% last modified: 22/03/2024
% 
% 20.07.01 update: adding logical for bathy > min model bathy threshold
% 20.07.31 updates:
% reading BATHY data before switch case doPlots YES/NO
% computing seaBedReview_temp out of switch case doPlots YES/NO
% 20.09.10 update: modifying computation of max depths measured by platform
% 21.01.18 updte: modifying computation of deepest measurement per profile
% 21.03.25 updates:
% modifying computation of bathymetry along platform trajectory
% (bathyTrip_temp) to avoid issues at edges of the map (additon of min/max
% in BATHY matrix index)
% modifying computation of deepest measurement per profile using TEMP
% instead of platform.TEMP (case of floats where resulting index
% mismatches)
% 22.02.03 update: add bathy in genData table
% 22.03.29 update: modifying computation of idx013 bathy min index (taking
% profile depth into consideration and extending to largest interval, case
% no location available)
% 22/03/2024 update: (jw)
% Script temporarily disabled: The ETOPO1 file (bathymetry shapefile) is too 
% large preventing pushing to GitHub, will aim to find a solution at some point
% (e.g. GIT LFS)
%
% -----------------------------------------------------------------------

disp(strcat('computing BATHYMETRY data along trajectory for platform:',...
    platform_metadata.platform_name,'.....'))


%% Read ETOPO file within the specified latLim and lonLim limits
% Extract bathymetry data within platform lat/lon limits from ETOPO dataset
ilatmin_in_bathy_temp = dsearchn(bathymetry.lat',floor(platform_metadata.geospatial_lat_min));
ilatmax_in_bathy_temp = dsearchn(bathymetry.lat',ceil(platform_metadata.geospatial_lat_max));
ilonmin_in_bathy_temp = dsearchn(bathymetry.lon',floor(convlon(platform_metadata.geospatial_lon_min,'signed')));
ilonmax_in_bathy_temp = dsearchn(bathymetry.lon',ceil(convlon(platform_metadata.geospatial_lon_max,'signed')));
if platform_metadata.geospatial_lon_max - platform_metadata.geospatial_lon_min < 180
    BATHY = bathymetry.data(ilatmin_in_bathy_temp:ilatmax_in_bathy_temp, ilonmin_in_bathy_temp:ilonmax_in_bathy_temp);
    BATHY_lon_temp = bathymetry.lon(ilonmin_in_bathy_temp : ilonmax_in_bathy_temp);
else
    % If the lon limits cross 180˚ E/W, the bathymetry subset needs to be pieced together
    BATHY_E_temp = bathymetry.data(ilatmin_in_bathy_temp:ilatmax_in_bathy_temp, ilonmax_in_bathy_temp:bathymetry.ref.RasterSize(2));
    BATHY_W_temp = bathymetry.data(ilatmin_in_bathy_temp:ilatmax_in_bathy_temp, 1:ilonmin_in_bathy_temp);
    BATHY = [BATHY_E_temp,BATHY_W_temp];
    BATHY_lon_E_temp = bathymetry.lon(ilonmax_in_bathy_temp:bathymetry.ref.RasterSize(2));
    BATHY_lon_W_temp = bathymetry.data(1:ilonmin_in_bathy_temp);
    BATHY_lon_temp = [BATHY_lon_E_temp,BATHY_lon_W_temp];
end
BATHY_lat_temp = bathymetry.lat(ilatmin_in_bathy_temp : ilatmax_in_bathy_temp);
lats_temp = BATHY_lat_temp';
lons_temp = BATHY_lon_temp';

%% PLOT WORLDMAP WITH FOCUS ON REGION OF DATASET + BATHYMETRY DATA
% TO BE REWRITTEN OR DELETED

switch doPlots
    case 'NO'
        
    case 'YES'
%         % Create figure
%         figure(131)
%         % Construct map axes for region
%         worldmap (platform_metadata.latLim,platform_metadata.lonLim)
% %         % alternatively (example)
% %         worldmap 'South America'
% 
% %         % Extract map projection structure from the current map axes
% %         mapProjStruct_temp = gcm ;
% %         mapAxes_temp = gca ;
% 
% %         % Define limits for lat and long according to worldmap region
% %         latLim = mapProjStruct_temp.maplatlimit ;
% %         lonLim = mapProjStruct_temp.maplonlimit ;
% 
% 
%         % Display map data, with extracted etopo value
%         h_temp = geoshow(BATHY, refvec_temp, 'DisplayType', 'texturemap') ;
% 
%         % Color the map based on the terrain elevation data, BATHY 
%         demcmap(BATHY, 500);
%         hc_temp = colorbar ;
%         s_temp = ['Elevation (m)'] ;
%         ylabel(hc_temp,s_temp) ;
% 
%         % Platform trip
%         plotm(platform.LATITUDE,platform.LONGITUDE,'.k','markersize',4)
% 
%         % Add marker on map (ex. for base station)
%         plotm(platform_metadata.base_station_LAT,...
%             platform_metadata.base_station_LON,...
%             'ok','markersize',4,'MarkerFaceColor','k')
%         % Naming on map
%         pointLabel_temp = platform_metadata.base_station ;
%         dx_temp = 0.1 ; dy_temp = 0.1 ; % displacement so the text does not overlay the data points
%         textm(platform_metadata.base_station_LAT + dx_temp,...
%             platform_metadata.base_station_LON + dy_temp,...
%             pointLabel_temp,...
%             'FontSize',8)
% 
% %         % set custom map limits
% %         setm(gca,'MapLatLimit',latLim)
% %         setm(gca,'MapLonLimit',lonLim)
% 
% %         % display coastlines
% %         load coastlines
% %         [latcells_temp, loncells_temp] = polysplit(coastlat, coastlon);
% %         numel(latcells_temp) ;
% %         plotm(coastlat, coastlon,'-k','LineWidth',0.5)
% %         % better use NOAA Shoreline / Coastline Resources (gshhg)
% 
%         % axes labels, figure title/legend
%         legend({strcat('platform:',platform_metadata.platform_name)},...
%             'Location','southeast')
%         hold off
%         % display figure on full screen
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]) ;
        
end


%% COMPUTE BATHYMETRY  + MAXIMUM MEASURED DEPTH FOR EACH PROFILE
% corresponding platform trajectory latitudes in BATHY LAT grid (lats_temp)
% (closest values selected)
idxlats_temp = dsearchn(lats_temp,platform.LATITUDE);
idxlats_temp(~idx001_latlonNonNan) = NaN ;

% corresponding platform trajectory longitudes in BATHY LON grid (lons_temp)
% (closest values selected)
idxlons_temp = dsearchn(lons_temp,platform.LONGITUDE);
idxlons_temp(~idx001_latlonNonNan) = NaN ;

% Compute bathymetry along platform trajectory
bathyTrip_temp = nan(np_tot,1) ;
for ii_temp = 1:np_tot
    if idx001_latlonNonNan(ii_temp) == 1
        bathyTrip_temp(ii_temp) = BATHY(idxlats_temp(ii_temp),idxlons_temp(ii_temp));
    end
end


% max depths measured by platform
idxDeepest_temp = find_ndim(~isnan(TEMP),1,'last').' ;
idxNull_temp = idxDeepest_temp == 0 ;
idxDeepest_temp(idxDeepest_temp == 0) = 1 ;
finiteIndices = find(isfinite(idxDeepest_temp));
% idxDeepest_temp = sub2ind(size(PRES),idxDeepest_temp,[1:platform_metadata.np].') ;
% deepestMeasPerProfile_temp = nan(platform_metadata.np,1) ;
% deepestMeasPerProfile_temp(~isnan(idxDeepest_temp)) =...
%     PRES(idxDeepest_temp(~isnan(idxDeepest_temp))) ;
% +++++++++++ FIXING SUB2IND ISSUE ++++++++++++++
deepestMeasPerProfile_temp = NaN(platform_metadata.np,1);
deepestMeasPerProfile_temp(finiteIndices,1) = - arrayfun(@(x) PRES(idxDeepest_temp(x), x), finiteIndices);
deepestMeasPerProfile_temp(idxNull_temp) = 0 ;
% deepestMeasPerProfile_temp = - deepestMeasPerProfile_temp ;

%% PLOT BATHYMETRY + MAX MEASURED DEPTHS

% highlight where measured platform depth > bathymetry database
seaBedReview_temp = deepestMeasPerProfile_temp ;
seaBedReview_temp(seaBedReview_temp >= bathyTrip_temp) = NaN ;
seaBedReview_temp(~idx001_latlonNonNan) = NaN ;


switch doPlots
    case 'NO'
        
    case 'YES'
        figure(132+1)
        plot(bathyTrip_temp,'-k')
        hold on
        plot(deepestMeasPerProfile_temp,'Color',[1 1 1] * 200/255) % light grey
        plot(seaBedReview_temp,'.r')
        hold off
        % edit axes
        xlabel('profile index')
        ylabel('depth (m)')
        title({'Vertical transect along platform trajectory',...
            strcat('platform:',platform_metadata.platform_name)})
        legend({...
            strcat('seabed elevation (',platform_metadata.bathyDataFileName(1:end),')'),...
            'platform max depth per profile',...
            'WARNING: platform max depth > bathymetry data'},...
            'Location','southoutside',...
            'Interpreter', 'none',...
            'Orientation','horizontal')
        % display figure on full screen
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]) ;

end



%% WRITE IN platform_metadata STRUCTURE

% bathymetry
platformBathy_temp = table(bathyTrip_temp,'VariableNames',{'bathy'}) ;

% profile max depth
platformBathy_temp.profDepth = deepestMeasPerProfile_temp ;

% logical for bathymetry
% (1) smaller than profile max depth i.e. to be revised
% (0) greater than profile max depth i.e. OK
seaBedReview_temp(~isnan(seaBedReview_temp)) = 1 ;
seaBedReview_temp(isnan(seaBedReview_temp)) = 0 ;
platformBathy_temp.seaBedReview = seaBedReview_temp ;

% logical for open ocean/coastal
% days in open ocean according to bathy data
daysOutBathy_temp =...
    genData.depDayNo(find(platformBathy_temp.bathy <...
    bathyMinModelThreshold)) ;
% days in open ocean according to profile depth
daysOutProfDepth_temp =...
    genData.depDayNo(find(platformBathy_temp.profDepth <...
    minProfDepthOpenOcean)) ;
% pick up largest interval
firstDayOut_temp = min(vertcat(daysOutBathy_temp,daysOutProfDepth_temp)) ;
lastDayOut_temp = max(vertcat(daysOutBathy_temp,daysOutProfDepth_temp)) ;
daysOut_temp = [firstDayOut_temp:lastDayOut_temp] ;
% create logical for open ocean=1/coastal=0
idxOpenOcean_temp = arrayfun(@(a) any(a == daysOut_temp),genData.depDayNo) ;
% diplay warning if empty vector
if nnz(idxOpenOcean_temp) == 0
    disp('*******************NO OPEN OCEAN PROFILES - REVIEW BATHY DATA*******************')
end
% write in platformBathy_temp table
platformBathy_temp.openOcean = idxOpenOcean_temp ;

platform_metadata.platform_bathy = platformBathy_temp ;


%%



%% create logical indexing for bathy >= bathyMinModelThreshold

idx013_bathyMin = platform_metadata.platform_bathy.openOcean ;


%% write bathy in genData table
genData.bathy = platform_metadata.platform_bathy.bathy ;


%% USE OF MAPPROFILE TO COMPUTE VALUES ALONG A TRANSECT

% mapprofile plots a profile of values between waypoints on a displayed regular data grid
% adapt from given example

% load korea
% [latlim, lonlim] = limitm(map, maplegend);
% figure
% worldmap(latlim, lonlim)
% meshm(map,maplegend,size(map),map)
% demcmap(map)
% plat = [40.5 30.7];
% plon = [121.5 133.5];
% [z,rng,lat,lon] = mapprofile(map,maplegend,plat,plon);
% plot3m(lat,lon,z,'w','LineWidth',2)
% figure
% plot(rng,z,'r')



%% CLEAR TEMP VARIABLES
clear *_temp


%% END
