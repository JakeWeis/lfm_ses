%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% load date array
% create numeric (dateNUM) and explicit (dateYMD) date arrays
% to be further used for scatter plots and solar position computation
% create genData table to store profile niformation
% and compute a few data about profile

%%%%%%%%%% REFERENCES %%%%%%%%%%
% function solarPosition
%   https://www.mathworks.com/matlabcentral/fileexchange/58405-solar-position-calculator
% function datepart
%   https://www.mathworks.com/matlabcentral/fileexchange/15585-datepart

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% J. Weis (jw)
% last modified: 22.03.29
% 
% 20.08.31 update: adusting DATE/PRES array size to max_depth_dataset value
% (in agreement with sesf012_gen_data_preprocess2_set_fixed_vert_grid
% script)
% 20.09.08 update: removing previous update for dateYMD array
% 20.09.23 updates:
% removing creation of PRES array (included in
% sesf012_gen_data_preprocess2_set_fixed_vert_grid)
% adding creation of genData table and computation of some data
% 22.02.02 updates:
% replacing genData.firstNonNan = NaN by genData.firstNonNan = 9999 case profile is all NaN
% adding test on profile depth/length in idx000_genProfileSelec
% adding a few data in genData table (date, lat/lon, platform name)
% 22.02.04 update: removing unused fields in genData table (MLDphySlm
% MLDt02 MLDphyBw
% 22.02.05 update: adding computation of day/night logical index (idx035)
% 
% -----------------------------------------------------------------------

disp(strcat('setting DATE data and computing basic PROFILE DATA for platform:',...
    platform_metadata.platform_name,'.....'))


%% load and set date arrays

% set reference date
dateref = '1950-01-01 00:00:00' ; % see platform.REFERENCE_DATE_TIME

% create dateYMD array for date in DD-MM-YYYY HH:MM:SS format
% (to be used for daylight filter)
dateYMD = datetime(platform.JULD_LOCATION + datenum(dateref),'ConvertFrom','datenum') ;
dateYMD = repmat(dateYMD',max_depth,1) ;

% create dateNUM array in numeric format
% (to be used for daylight filter)
dateNUM = platform.JULD_LOCATION + datenum(dateref) ;
dateNUM = dateNUM.' ; % transpose matrix to get a 1 x n array (display homogeneity)


%% create genData table to store general profile information
% short names
names_temp = {} ;
names_temp{1} = 'profileNB' ; % profile number
names_temp{2} = 'firstNonNan' ; % depth of first non NaN value of light
names_temp{3} = 'lastNonNan' ; % depth of last non NaN value of light
names_temp{4} = 'hOfDayUTC' ; % hour of the day (UTC) of the profile
names_temp{5} = 'solAlt' ; % solar altitude value (angle, in dregrees)
names_temp{6} = 'MLDphy' ; % value to be used for MLDphy (either MLDphyDHT2009 or MLD003)
names_temp{7} = 'MLDHTd2009' ; % MLD from Holte and Talley 2009 algorithm
names_temp{8} = 'MLD003' ; % MLD based on threshold criteria deltad = 0.03 kg.m-3


genData = horzcat([1:np_tot].', nan(np_tot,numel(names_temp)-1)) ;
genData = array2table(genData,'VariableNames',{names_temp{1:end}}) ;


%% Compute a few data about profiles and write in genData table

% LAT/LON
genData.LONGITUDE = platform.LONGITUDE ;
genData.LATITUDE = platform.LATITUDE ;

% date NUM/YMD / deployment day Nº
genData.dateNUM = dateNUM.' ;
genData.dateYMD = dateYMD(1,:).' ;
genData.depDayNo = (floor(dateNUM) - floor(dateNUM(1)) + 1).' ;

% platform name
genData.platform_name = repmat(platform_metadata.platform_name,platform_metadata.np,1) ;

% first/last nonNaN value of each profile (based on TEMP field)
genData.firstNonNan = find_ndim(~isnan(platform.TEMP),1,'first').' ;
genData.firstNonNan(genData.firstNonNan == 0) = 9999 ;
genData.lastNonNan = find_ndim(~isnan(platform.TEMP),1,'last').' ;
genData.lastNonNan(genData.lastNonNan == 0) = 9999 ;

% write false in idx000_genProfileSelec where profiles are all NaN
idx000_genProfileSelec(genData.firstNonNan > max_depth) = false ;

% write false in idx000_genProfileSelec where profiles do not meet profile length criteria
idx000_genProfileSelec(genData.lastNonNan <= minProfileDepth) = false ;
idx000_genProfileSelec(genData.lastNonNan - genData.firstNonNan < minProfileLength) = false ;

% parameters for computation of solar altitude
rot_temp = 0;       % [arc-degrees] rotation clockwise from north
zen_az_temp = [0 0] ;

for ii_temp = 1:np_tot
%     if isnan(genData.firstNonNan(ii))
%     else  
        % compute solar angle and solar altitude (even for all NaN profiles)   
        % solar angle
        if idx001_latlonNonNan(ii_temp) == 0
            zen_az_temp = [NaN NaN] ;
        else
            zen_az_temp = solarPosition(dateNUM(ii_temp),...
                platform.LATITUDE(ii_temp),platform.LONGITUDE(ii_temp),...
                TimOffZ,rot_temp,DaySavTim) ;
        end    
        % hour of the day (UTC)
        Hdec_temp = datepart(dateYMD(1,ii_temp),'hour') ;
%     end
	% write computed values in genData array
    genData.solAlt(ii_temp) = 90 - zen_az_temp(1) ;	% solar altitude value (angle, in dregrees)
	genData.hOfDayUTC(ii_temp) = Hdec_temp ;         % hour of the day (UTC) of the profile   
end


%% create logical indexing for daylight threshold = xxº solar altitude

% daylight profiles
idx035_lightDay = genData.solAlt >= sol_elv_threshold ; % logical for solar altitude >= 20º (daylight filter)
disp(strcat('DAYLIGHT THRESHOLD =',32,num2str(sol_elv_threshold),'º (Solar Altitude)'))


%% CLEAR TEMP VARIABLES
clear *_temp


%% END
