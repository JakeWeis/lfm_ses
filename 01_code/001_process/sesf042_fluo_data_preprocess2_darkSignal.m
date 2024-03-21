%% script metadata
% alternative script with merged per-profile / per-platform approach

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% estimate FLUO dark signal (profile by profile / for each platform)
% correct FLUO with dark signal (profile by profile)
% remark: in the following script, the "dark" signal is the part of the signal
% where FLUO is supposed to be zero i.e. at depth (maxCHLA_depth)

%%%%%%%%%% REFERENCES %%%%%%%%%%
% see Bellacicco et al. 2019 doi: 10.1029/2019GL084078 (Supplementary
% Information)

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x
% sesf01x
% sesf04x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% J. Weis (jw)
% last modified: 22.03.29

% 20.02.10 update:
% DARK signal = cst part of profile
% no correction applied to rest of profile (DARK signal set to NaN)
% 20.03.31 update:
% add histogram of chla dark values
% write DARK value in chlaData_array
% 20.03.31 update: adding switch case for plot/no plot option (default =
% YES)
% 20.08.28 update: adding condition in histogram plot case CHLA dataset is all NaN
% 21.01.12 update: changing _adReg names to _nadReg (use of non ADJUSTED
% CHLA)
% 21.01.18 update:
% including replacement by NaN of Chl-a data for depths > maxCHLA_depth
% corresponding values in chlaDatatable are majored by maxCHLA_depth
% 21.04.26 update: updating computation of dark signal based on merged
% per-profile / per-platform approach
% 22.03.29 update: chla becomes fluo (names of variables) e.g. creation of
% FLUO_nadRegDk (previously named CHLA_nadRegDk)
% 
% -----------------------------------------------------------------------

disp(strcat('FLUO data pre-processing step2: dark signal / platform:',...
    platform_metadata.platform_name,'.....'))

%% compute dark value on the last delta meters of the profile (default = 10)

% interval to compute median value of deep measurements
darkDepthDelta_temp = 10 ;
% array of interest
fluoDeep_temp = FLUO_nadReg(maxCHLA_depth - (darkDepthDelta_temp - 1) :...
    maxCHLA_depth,:) ;
fluoDeep_temp(:,~idx041_chlaNonNan) = NaN ;
% median value of deep measurements per profile (if profile reaches depth)
fluoProfDarkVal_temp = median(fluoDeep_temp,1) ;
% fill missing values with closest non NaN value
fluoProfDarkVal_temp = fillmissing(fluoProfDarkVal_temp,'nearest') ;

% write dark value in chlaData table
fluoData.darkValProfile = fluoProfDarkVal_temp.' ;

% smooth dark value (to be used for correction)
% time window to be used for smoothing (days)
nDays_temp = 10 ;
% apply moving median (not computed in edges of dataset)
darkValCorr_temp = movmedian(fluoData.darkValProfile,...
    platform_metadata.np / (1/nDays_temp * max(genData.depDayNo)),...
    'Endpoints','fill') ;
% fill missing values with closest value
darkValCorr_temp = fillmissing(darkValCorr_temp,'nearest') ;
% fit dark values to approx. by linear drift (as function of time)
darkValCorrFit_temp = fitlm(1:platform_metadata.np,darkValCorr_temp) ;
fluoData.darkValCorr = darkValCorrFit_temp.Fitted ;

%% remove dark value offset for each profile

% initialize CHLA_nadRegDk array
FLUO_nadRegDk = FLUO_nadReg ;
% offset to be applied to CHLA measurement
fluoOffset_temp = repmat(fluoData.darkValCorr,1,max_depth).' ;
% remove offset
FLUO_nadRegDk = FLUO_nadRegDk - fluoOffset_temp ;
% erase negative values
FLUO_nadRegDk(FLUO_nadRegDk < 0) = 0 ;


%% detect starting depth of cst FLUO signal
% (if occurs due to e.g. electronic/data processing artefact)
% + remove corresponding part of the profile

for ii_temp = find(idx041_chlaNonNan).'
    % initialize PDARK_temp at each loop
    PDARK_temp = maxCHLA_depth + 1 ;
    % index of first/last non NaN value of FLUO
    firstNonNan_temp = fluoData.firstNonNan(ii_temp) ;
    lastNonNan_temp = fluoData.lastNonNan(ii_temp) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % detect FLUO cst values (constant tail of profile)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Sampled vector X must have at least 2 valid observations
    if lastNonNan_temp - firstNonNan_temp + 1 < 2
        PDARK_temp = NaN ;
    else
        for jj_temp = firstNonNan_temp:lastNonNan_temp - 1
            % proceed to cst test (dark offset has been previously removed)
            amplitude_temp = nanmax(FLUO_nadRegDk(jj_temp:lastNonNan_temp,ii_temp)) ;
            if amplitude_temp <= delta_for_cstCHLA
                break
                % end for loop when range < delta_for_cstCHLA i.e. FLUO signal is considered constant in interval j:botom
            end
        end      
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % correct FLUO profile (i.e. erase cst values)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if jj_temp == lastNonNan_temp - 1
            % no values to discard if amplitude_temp never proves to be <
            % delta_for_cstCHLA
        else
            % discard cst values
            ind_dark_temp = jj_temp ;
            PDARK_temp = PRES(ind_dark_temp,ii_temp) ;
            % set cst part to zero
            FLUO_nadRegDk(ind_dark_temp:lastNonNan_temp,ii_temp) = 0 ;
        end
    end
    % write computed values in chlaData array
    fluoData.darkDepth(ii_temp) = PDARK_temp ;

end

%% plot dark value data

switch doPlots
    case 'NO'
        
    case 'YES'

        % plot time series of FLUO dark values with interval length and histogram
        % of DARK values
        figure(2201)
        sgtitle(strcat('time series of FLUO dark value - tag:',platform_metadata.platform_name))

        subplot(2,1,1)
        plot([1:np041_chlaNonNan],...
            fluoData.darkValProfile(idx041_chlaNonNan),'.k',...
            'DisplayName','dark value (per profile)')
        hold on
        plot([1:np041_chlaNonNan],...
            fluoData.darkValCorr(idx041_chlaNonNan),'-k',...
            'DisplayName','effective dark correction (movmedian + linear fit)')
        xlabel('profile nº')
        ylabel('FLUO dark value (mg.m-3)')
        hc_temp = colorbar ;
        hc_temp.Location = 'eastoutside';
        ylabel(hc_temp,'profile interval length (m)');
        legend('Location','northwest')

        subplot(2,1,2)
        histogram(fluoData.darkValProfile(idx041_chlaNonNan),...
            'DisplayName','histogram of FLUO dark value')
        X_temp = fluoData.darkValProfile(idx041_chlaNonNan) ;
        p_temp = 95 ; % eventually adjust percentile to be shown on plot
        Y_temp = quantile(X_temp,p_temp / 100) ;
        % xlim([0 floor(3*Y_temp)])
        xlabel('FLUO dark value (mg.m-3)')
        ylabel('frequency')
        if ~isnan(Y_temp)
            h_temp = xline(Y_temp,...
                'DisplayName',strcat(num2str(p_temp,2),'th percentile')) ;
            h_temp.LineStyle = '--' ;
        end
        legend('Location','NorthEast')

        % set figure size to whole screen height
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 1]) ; 

end


%% clear temp_data

clear *_temp

%% END
