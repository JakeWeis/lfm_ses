%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% preprocess PAR data
% detect PAR dark signal (PAR "absolute zero")
% correct PAR profile with dark value (profile by profile)

%%%%%%%%%% REFERENCES %%%%%%%%%%
% Organelli et al. 2016

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x
% sesf01x
% sesf03x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% last modified: 22.04.06
% 
% 20.01.30 update: no more correction, just nan after dark value detected
% 20.02.10 update: insert removing of constant tail of profile before normality test
% 20.02.10 update: (linear) interpolation after removal of negative values due to dark correction
% (keep positivity, avoid Epsilon setting issues)
% 20.04.07 update: adding calculation of lightData_ad.attSlopPART
% 20.04.21 update: modifying filling with nan in dk corrected vector
% 20.04.27 update: fixing issues with fillmissing function (re-introducnig negative values)
% 20.09.10 updates: small adjustments on indexes to be erased by nan in
% constant part of the profile
% 21.01.12 update: removing if test on idx031_lightNonNan
% 21.01.18 update: remove displaying of processing step
% 21.04.27 update:
% changing LIGHT array variable name
% L_adNonLogRegDkCted into L_adNonLogRegDk for naming consistency
% adding computation of L_adLogRegDk
% 22.02.02 updates:
% adding 'EndValues','none' option to fillmissing when
% replacing negative values by NaN (end of dark correction process) to
% avoid re-writing of negative values (leading to creating complex numbers
% when computing log)
% compute dark depth only for the 5 deepest profiles per deployment day
% 22.02.03 update: separate dark noise detection (profile per profile / n
% deepest profiles per day) and dark noise correction (matricial)
% 22.04.06 update: rename LIGHT into PAR (L_adxxx becomes PAR_nadxxx)
%
% -----------------------------------------------------------------------

disp(strcat('PAR data pre-processing step2: dark signal / platform:',...
    platform_metadata.platform_name,'.....'))

%% initialize variables

% create PAR_nadNonLogDkCted array corresponding to PAR_nadNonLog corrected by calculated dark value
PAR_nadNonLogRegDk = PAR_nadNonLogReg ;
% PAR_nadNonLogRegDk_temp computed to display dark noise - temp data
PAR_nadNonLogRegDk_temp = PAR_nadNonLogReg ;


%% define subset of processed profiles

% select a few profiles per day for dark correction: the nProfiles deepest per day
% nProfilesPerDay_dark is set in script sesf004

% starting index of each deployment day Nº
[unVal_temp,unInd_temp] = unique(genData.depDayNo) ;
% compute 5 deepest per day
[deepestPerDayVal_temp,deepestPerDayInd_temp] =...
    arrayfun(@(a) maxk(parData.lastNonNan(genData.depDayNo == a),nProfilesPerDay_dark),...
    unVal_temp,'UniformOutput',false) ;
% correction of the index to take into account 1st index of each deployment day
deepestPerDayIndCorr_temp = arrayfun(@(a) deepestPerDayInd_temp{a} + unInd_temp(a) - 1,...
    1:numel(unVal_temp),'UniformOutput',false) ;
% concatenate indexes in one array
deepestPerDayIndCorr_temp = horzcat(deepestPerDayIndCorr_temp{1:end}) ;
deepestPerDayIndCorr_temp = deepestPerDayIndCorr_temp(:) ;
% sort indexes
deepestPerDayIndCorr_temp = sort(deepestPerDayIndCorr_temp,'ascend') ;

% resulting computation index
indexDeepest_temp = false(platform_metadata.np,1) ;
indexDeepest_temp(deepestPerDayIndCorr_temp) = true ;
indexes_temp = indexDeepest_temp & idx031_lightNonNan ;

% display computation progress every 25 analyzed profiles
dispInfo_temp = find(indexes_temp) ;
dispInfo_temp = dispInfo_temp(1:25:end) ;

%% let's go: detect PAR zero dark value (normality test, organelli 2016)

tic % start stopwatch timer


for iProfile_temp = find(indexes_temp).'
    % display computation progress every 100
    if nnz(iProfile_temp == dispInfo_temp) == 1
        disp(strcat('ddepth platform:',32,char(platform_metadata.platform_name),...
            '/profile:',32,int2str(iProfile_temp),'/',int2str(platform_metadata.np)))
    end
    % initialize PDARK_temp and VDARK_temp at each loop
    PCST_temp = max_depth ;
    PDARK_temp = max_depth ;
    VDARK_temp = 0 ;
    first_non_NaN_temp = parData.firstNonNan(iProfile_temp) ;
    last_non_NaN_temp = parData.lastNonNan (iProfile_temp) ;
        
        
        
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % detect (eventual) constant part of the signal -> data range = Espilon
    % (remove tail of profile to avoid normality test on constant data series - not working well)
    % example:
    % [h_temp,p_temp] = lillietest(ones(1,max_depth),'Alpha',0.01,'Distribution','normal') ;
    % returns h_temp = 1 (i.e. rejection of the null hypothesis)
    % null hypothesis = distribution is normal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%


    % initialize cstDetect_temp used in cst test
    cstDetect_temp = 0 ;        
    % Sampled vector X must have at least 2 valid observations
    if last_non_NaN_temp - first_non_NaN_temp + 1 < 2
        ind_cst_temp = last_non_NaN_temp + 1 ;
        PCST_temp = NaN ;
    else
        for jj_temp = first_non_NaN_temp:last_non_NaN_temp - 1 % Sampled vector X must have at least 2 valid observations
            % proceed to cst test
            range_temp = range(PAR_nadNonLogRegDk(jj_temp:last_non_NaN_temp,iProfile_temp)) ;
            if range_temp <= delta_for_cstLIGHT
                cstDetect_temp = 1 ;
                break
                % end for loop when range < delta_for_cstLIGHT i.e. PAR signal is considered constant in interval j:botom
            end
        end

        % write index of first value belonging to dark signal
        ind_cst_temp = jj_temp + 1 ;
        % write depth of dark value
        PCST_temp = PRES(ind_cst_temp,iProfile_temp) ;

    end

    % set cst part to NaN
    switch cstDetect_temp
        case 0
            % do not affect PAR_nadNonLogRegDkCted(:,i) vector
        case 1
            % (use of min to exclude case where ind_cst_temp corresponds to
            % last value of profile)
            PAR_nadNonLogRegDk(min(ind_cst_temp + 1,max_depth):end,iProfile_temp) = NaN ;
    end





    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % detect PAR zero dark value (normality test, organelli 2016)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%



    % set new last non NaN value after detection of cst tail of signal
    switch cstDetect_temp
        case 0
            % case of no cst detection: keep last value of profile
        case 1
            % case of cst detection: take first value belonging to dark
            % signal - 1
            last_non_NaN_temp = ind_cst_temp - 1 ;
    end
    % initialize h_temp used in lillietest
    h_temp = 1 ;
    % Sampled vector X in Lilliefors test must have at least 4 valid observations
    if last_non_NaN_temp - first_non_NaN_temp + 1 < 4
        ind_dark_temp = last_non_NaN_temp + 1 ;
        PDARK_temp = NaN ;
        VDARK_temp = NaN ;
        slop_temp = NaN ;
        slopPART_temp = NaN ;
    else
        for jj_temp = first_non_NaN_temp:last_non_NaN_temp - 3 % Sampled vector X must have at least 4 valid observations
            % turn off following warning message
            % Warning: P is less than the smallest tabulated value, returning 0.001.
            warning('off','stats:lillietest:OutOfRangePLow')
            warning('off','stats:lillietest:OutOfRangePHigh')
            % proceed to Lilliefors test
            [h_temp,p_temp] = lillietest(PAR_nadNonLogRegDk(jj_temp:last_non_NaN_temp,iProfile_temp),'Alpha',0.01,'Distribution','normal') ;
            % use name-value pair argument 'MCTol',1e-4 to find a more accurate p-value

            if h_temp == 0
                break
                % end for loop when h=0 i.e. distribution of PAR values in interval j:botom is normal
            end
        end

        % write index of first value belonging to dark signal
        switch h_temp
            case 1 % case of all null lillietest: keep last value of profile
                ind_dark_temp = last_non_NaN_temp ;
                % keep depth of dark value at last non nan measurement
                % of profile and dark value at zero (default)
                PDARK_temp = PRES(ind_dark_temp,iProfile_temp) ;
%                     VDARK_temp = 0 ;

            case 0 % case of lillietest detection: keep index where loop broke
                ind_dark_temp = jj_temp ;
                % write dark value and depth of dark value
                PDARK_temp = PRES(ind_dark_temp,iProfile_temp) ;
                VDARK_temp = median(PAR_nadNonLogRegDk(ind_dark_temp:last_non_NaN_temp,iProfile_temp)) ;
        end


        % slope of attenuation (calculated where PAR is non NaN / non DARK)
        slop_temp = (PAR_nadLogReg(ind_dark_temp,iProfile_temp) -...
            PAR_nadLogReg(first_non_NaN_temp,iProfile_temp)) /...
            (PRES(ind_dark_temp,iProfile_temp) -...
            PRES(first_non_NaN_temp,iProfile_temp)) ;
        % slope of attenuation
        % (calculated between min_topIntegBound and min_botIntegBound)
        if ind_dark_temp > min_botIntegBound
            slopPART_temp = (PAR_nadLogReg(min_botIntegBound,iProfile_temp) -...
                PAR_nadLogReg(first_non_NaN_temp,iProfile_temp)) /...
                (PRES(min_botIntegBound,iProfile_temp) -...
                PRES(first_non_NaN_temp,iProfile_temp)) ;
        else
            slopPART_temp = NaN ;
        end

    end

    % write computed values in parData array
    parData.attSlopTOT(iProfile_temp)  = slop_temp ;	% slope of attenuation (calculated where PAR is non DARK)
    parData.attSlopPART(iProfile_temp) = slopPART_temp ; % slope of attenuation calculated on [min_topIntegBound:min_botIntegBound] interval
    parData.darkDepth(iProfile_temp)   = PDARK_temp ;	% starting depth of dark value (organelli 2016 method)
    parData.darkVal(iProfile_temp)     = VDARK_temp ;	% dark value (organelli 2016 method)

end   
        
        
%% extend computed values to all profiles

% remove zeros
parData.darkDepth(parData.darkDepth == 0) = NaN ;
parData.darkVal(parData.darkVal == 0) = NaN ;

% fill missing values
parData.darkDepth = fillmissing(parData.darkDepth,'nearest') ;
parData.darkVal = fillmissing(parData.darkVal,'nearest') ;


%% correct PAR profiles with dark value

% compute mask of dark signal
darkNaNMask_temp = arrayfun(@(a) vertcat(false(a - 1,1),true(max_depth - a + 1,1)),...
    parData.darkDepth.','UniformOutput',false) ;
darkNaNMask_temp = horzcat(darkNaNMask_temp{1:end}) ;

% set dark noise to NaN according to dark mask
PAR_nadNonLogRegDk(darkNaNMask_temp) = NaN ;

% compute mask of dark values
darkValMask_temp = repmat(parData.darkVal.',max_depth,1) ;

% correct PAR_nadNonLogRegDkCted by calculated dark values
PAR_nadNonLogRegDk = PAR_nadNonLogRegDk - darkValMask_temp ;

% remove remaining negative values in effective signal
negidx_temp = find(PAR_nadNonLogRegDk <= 0) ;
PAR_nadNonLogRegDk(negidx_temp) = NaN ; % replace negative values by NaN

% replace nan in vector profile using linear interpolation
% replace only if more than 1 non nan value in raw dataset
% number of non nan values in vect_to_reg_temp
n_values_temp = arrayfun(@(a) nnz(~isnan(PAR_nadNonLogRegDk(:,a))),...
    [1:platform_metadata.np],'UniformOutput',false) ;
n_values_temp = horzcat(n_values_temp{1:end}).' ;
PAR_nadNonLogRegDk(:,n_values_temp >= 2) = fillmissing(PAR_nadNonLogRegDk(:,n_values_temp >= 2),...
    'linear',1,'EndValues','none') ;


%% PAR_nadNonLogDkCted_temp is computed to display dark noise if necessary
% Nan are replaced by zeros
PAR_nadNonLogRegDk_temp = PAR_nadNonLogRegDk_temp - darkValMask_temp ;
PAR_nadNonLogRegDk_temp(negidx_temp) = 0 ;



toc % stop stopwatch timer

%% end of timer: breathe a bit



%% compute mean of dark value (excluding NaN values) - to be compared to Epsilon
mean_VDARK_temp = nanmean(parData.darkVal(:))


%% compute natural log of dark-corrected PAR data
PAR_nadLogRegDk = log(PAR_nadNonLogRegDk) ;


%% clear temp data

clear ans *_temp

disp('script with dark correction performed')

%% END
    