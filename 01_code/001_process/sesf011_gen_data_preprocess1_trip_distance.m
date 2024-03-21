%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% calculate cumulated distance from beginning of trip and relative distance
% for each profile
% when LAT/LON is available
% store in platform_metadata structure

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x
% sesf01x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
disp('***MATLAB version information***')
disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% last modified: 22.03.29

% 20.02.12 update: add colorbar on trajectory plot to display date
% 20.03.31 update: add histogram of relative distance
% 20.03.31 update: adding switch case for plot/no plot option (default =
% YES)
% 22.03.29 update: adding grid on tajectory plot
% 
% -----------------------------------------------------------------------

disp(strcat('processing TRIP distance data for platform:',...
    platform_metadata.platform_name,'.....'))


%% display trajectory of platform

switch doPlots
    case 'NO'
        
    case 'YES'
        figure(111)
        scatter(platform.LONGITUDE(idx001_latlonNonNan),...
            platform.LATITUDE(idx001_latlonNonNan),...
            [],dateNUM(idx001_latlonNonNan),'.')
        % add base stationn location
        hold on
        plot(platform_metadata.base_station_LON,platform_metadata.base_station_LAT,'ok')
        dx_temp = 0.01 ; dy_temp = 0.01; % displacement so the text does not overlay the data points
        text(platform_metadata.base_station_LON + dx_temp,...
            platform_metadata.base_station_LAT + dy_temp,...
            platform_metadata.base_station);
        hold off
        % add labels, colorbar and edit plot properties
        ylabel('LAT')
        % lat_min = floor(min(platform.LATITUDE)) ;
        % lat_max = floor(max(platform.LATITUDE)) + 1 ;
        % ylim([lat_min lat_max])
        % yticks([lat_min:1:-lat_max])
        xlabel('LON')
        % lon_min = floor(min(platform.LONGITUDE)) ;
        % lon_max = floor(max(platform.LONGITUDE)) + 1 ;
        % xlim([lon_min lon_max])
        % xticks([lon_min:1:lon_max])
        axis equal
        title({'this was my trip',strcat('platform:',platform_metadata.platform_name)})
        hc_temp = colorbar ;
        hc_temp.Location = 'eastoutside' ;
        s_temp = ['date'] ;
        ylabel(hc_temp,s_temp) ;
        date_ticks_temp = get(hc_temp,'XTick') ;
        date_ticks_temp = datestr(date_ticks_temp,'dd/mm/yy') ;
        set(hc_temp,'XTickLabel',date_ticks_temp)
        colormap(gca,'parula')
        grid on

end

%% calculate cumulated/relative distance at date D of profile i

% create platform_distancekm array for cumulated/relative distance computation
% such that

% platform_distancekm array
% col 1 = column k is the cumulated distance at date platform.JULD_LOCATION(k)
% col 2 = column k is the relative distance to previous profile(k-1)

platform_distancekm_temp = zeros(np_tot,2);

for jj_temp = 2:np_tot
    if and(idx001_latlonNonNan(jj_temp-1),idx001_latlonNonNan(jj_temp)) == 0
        platform_distancekm_temp(jj_temp,1) = platform_distancekm_temp(jj_temp-1,1) ;
        platform_distancekm_temp(jj_temp,2) = NaN ;
    else
        point1_temp = [platform.LATITUDE(jj_temp-1) platform.LONGITUDE(jj_temp-1)] ;
        point2_temp = [platform.LATITUDE(jj_temp) platform.LONGITUDE(jj_temp)] ;
        [dist1km_temp dist2km_temp] = lldistkm(point1_temp, point2_temp) ;
        % write computed values in platform_distancekm array
        platform_distancekm_temp(jj_temp,1) = platform_distancekm_temp(jj_temp-1,1) + dist1km_temp ;  % col 1 = cumulated distance
        platform_distancekm_temp(jj_temp,2) = dist1km_temp ;	% col 2 = relative distance
    end
end

% convert to table with line names
platform_distancekm_temp =...
    array2table(platform_distancekm_temp,...
    'VariableNames',{'CUM','REL'}) ;

platform_metadata.platform_distance = platform_distancekm_temp ;


%% plot cumulated distance and relative distance

switch doPlots
    case 'NO'
        
    case 'YES'
        figure (112)
        sgtitle(strcat('platform:',platform_metadata.platform_name,...
                        ' (', int2str(np_tot),' profiles)'))

        % cumulated distance
        subplot(3,1,1)
        plot(dateYMD(1,idx001_latlonNonNan),platform_metadata.platform_distance.CUM(idx001_latlonNonNan))
        xlabel('date')
        ylabel('cumulated distance (km)')
        legend(strcat('TOTAL=',int2str(floor(platform_metadata.platform_distance.CUM(end))),' km'),'Location','SouthEast')

        % relative distance
        subplot(3,1,2)
        plot(dateYMD(1,idx001_latlonNonNan),platform_metadata.platform_distance.REL(idx001_latlonNonNan))
        xlabel('date')
        ylabel('relative distance (km)')

        % histogram
        subplot(3,1,3)
        histogram(platform_metadata.platform_distance.REL(idx001_latlonNonNan))
        X_temp = platform_metadata.platform_distance.REL(idx001_latlonNonNan) ;
        p_temp = 99 ; % eventually adjust percentile to be shown on plot
        Y_temp = quantile(X_temp,p_temp / 100) ;
        xlim([0 floor(3*Y_temp)])
        xlabel('relative distance (km)')
        ylabel('frequency')
        h_temp = xline(Y_temp) ;
        h_temp.LineStyle = '--' ;
        legend({'histogram of relative distance',...
            strcat(num2str(p_temp,2),'th percentile')},...
            'Location','NorthEast')

        % set figure size to whole screen height
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 1]) ; 

end


%% clear temp data
clear *_temp

%% END