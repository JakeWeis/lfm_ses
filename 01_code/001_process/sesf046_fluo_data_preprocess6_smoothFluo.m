%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% smooth FLUO profile
% METHOD1
% moving median with 5 m window followed by 7-point mean filter
% (from Lacour et al. 2017 additional data)
% 
% METHOD2
% proceed to functional fit
% 2 fits are computed:
%     - 1 fit on a fixed interval (min_topIntegBound:min_botIntegBound)
%     - 1 fit on the entire vertical profile where data is non Nan

%%%%%%%%%% REFERENCES %%%%%%%%%%
% fda package
% https://www.psych.mcgill.ca/misc/fda/software.html

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x
% sesf01x
% sesf04x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
disp('***MATLAB version information***')
disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% last modified: 22.03.29

% 20.03.31 update: adding switch case for plot/no plot option (default =
% YES)
% 20.11.19 update: all plots on same plot (hold on) for visual check
% (enabled if doPlots = 'YES')
% 21.01.12 update:
% adjusting variables names (CHLA_nad and chlaData)
% adding condition in lasNonNan index with max value = maxCHLA_depth
% 21.04.27 updates:
% replacing Butterworth filter by moving median + moving mean filter
% performing smoothing after NPQ correction (done before NPQ in previous
% version)
% 21.07.21 update: replacing fit by functional fit (fda)
% 
% -----------------------------------------------------------------------

disp(strcat('FLUO data pre-processing step6: Smoothing / platform:',...
    platform_metadata.platform_name,'.....'))


%%%%% METHOD1 %%%%%
%% smooth Chl-a data - moving median / moving mean

% % 5-point moving median
% fluoSthMed5_temp = movmedian(FLUO_nadRegDkNpq,5,1) ;
% 
% % 7-point moving mean
% fluoSthMean7_temp = movmean(fluoSthMed5_temp,7,1) ;


%%%%% METHOD2 %%%%%
%% initialize FLUO_nadRegDkNpqSth array

% fitted data (entire interval)
fluoFdaFitAll_temp = nan(max_depth,platform_metadata.np) ;
FLUO_nadRegDkNpqFitAll = nan(max_depth,platform_metadata.np) ;

% fitted data (fixed interval)
fluoFdaFitBnd_temp = nan(max_depth,platform_metadata.np) ;
FLUO_nadRegDkNpqFitBnd = nan(max_depth,platform_metadata.np) ;

% fd objects
fluo_fdAllCoefs = nan(nbaz,platform_metadata.np) ;
fluo_fdBndCoefs = nan(nbaz,platform_metadata.np) ;

%% COMPUTE FIT ON ENTIRE VERTICAL PROFILE
% SELECT PROFILES TO COMPUTE
fluoToFitBnd_temp = FLUO_nadRegDkNpq ;
indexes_temp =...
    idx000_genProfileSelec &...
    idx041_chlaNonNan ;  
indFitAll_temp = find(indexes_temp) ;

% DATA2FD PARAMETERS
% define basis
nbaz_temp = nbaz ;
nord_temp = nord ;

% FIT FLUO
% compute fit
for iProfile_temp = indFitAll_temp.'
    % parameters
    vec0_temp = fluoToFitBnd_temp(:,iProfile_temp) ;
    pp_temp = find(~isnan(vec0_temp)) ;
    vec_temp = vec0_temp(pp_temp) ;
    wbasis_temp = create_bspline_basis([pp_temp(1) pp_temp(end)],nbaz_temp,nord_temp) ;
    % fit
    vecfd0_temp = data2fd(pp_temp,vec_temp,wbasis_temp) ;
    % eval fit
    evalvecfd_temp = eval_fd(vecfd0_temp,pp_temp) ;
    fluoFdaFitAll_temp(pp_temp,iProfile_temp) = evalvecfd_temp ;
    fluo_fdAllCoefs(:,iProfile_temp) = getcoef(vecfd0_temp) ;
end


%% COMPUTE FIT ON A FIXED INTERVAL
% SELECT PROFILES TO COMPUTE
fluoToFitBnd_temp = FLUO_nadRegDkNpq ;
idxBnd_temp =...
    find_ndim(~isnan(fluoToFitBnd_temp),1,'first') <= min_topIntegBound &...
    find_ndim(~isnan(fluoToFitBnd_temp),1,'last') >= min_botIntegBound ;
idx046_chlaFDAfitBnd =...
    idxBnd_temp.' &...
    idx000_genProfileSelec &...
    idx041_chlaNonNan ;
indFitBnd_temp = find(idx046_chlaFDAfitBnd) ;

% DATA2FD PARAMETERS
% depth interval
pp_temp = (min_topIntegBound:1:min_botIntegBound).' ;
% Data to fit/convert to fd object
fluoToFd_temp = fluoToFitBnd_temp(pp_temp,indFitBnd_temp) ;
% define basis
nbaz_temp = nbaz ;
nord_temp = nord ;
wbasis_temp = mybL ;

% FIT
vecfdBnd_temp = data2fd(pp_temp,fluoToFd_temp,wbasis_temp) ;

% EVAL FIT
fluoFdaFitBnd_temp(pp_temp,indFitBnd_temp) = eval_fd(vecfdBnd_temp,pp_temp) ;
fluo_fdBndCoefs(:,indFitBnd_temp) = getcoef(vecfdBnd_temp) ;


%% write data
% FLUO_nadRegDkNpqSth = fluoSthMean7_temp ;
FLUO_nadRegDkNpqFitAll = fluoFdaFitAll_temp ;
FLUO_nadRegDkNpqFitBnd = fluoFdaFitBnd_temp ;


%% visual check

switch doPlots
    case 'NO'
        
    case 'YES'
        figure (2301)
        % full screen
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]) ;
        for ii_temp = find(idx041_chlaNonNan).'
            plot(FLUO_nadReg(:,ii_temp),PRES(:,ii_temp),'.k',...
                'DisplayName','raw')
            hold on
            plot(FLUO_nadRegDkNpq(:,ii_temp),PRES(:,ii_temp),'og',...
                'DisplayName','dark')
            plot(FLUO_nadRegDkNpq(:,ii_temp),PRES(:,ii_temp),'-g',...
                'DisplayName','dark+NPQ')
            plot(FLUO_nadRegDkNpqFitAll(:,ii_temp),PRES(:,ii_temp),'-k',...
                'DisplayName','dark+NPQ+smoothed')
            set(gca,'YDir','reverse')
            xlabel('FLUO (mg.m-3)')
            ylabel('depth (m)')
            title(strcat('Chlorophyll a - platform:',platform_metadata.platform_name),...
                strcat(' - profile nº',...
                int2str(nnz(idx041_chlaNonNan(1:ii_temp))),...
                ' /',int2str(np041_chlaNonNan)))
            legend('Location','eastoutside')
            pause(0.5)
            hold off
        end
end


%% clear temp_data

clear *_temp


%% END