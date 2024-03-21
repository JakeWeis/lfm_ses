%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% load FLUO/LIGHT/TEMP/PSAL arrays
% fit dataset to a fixed vertical grid (1 line <-> 1 meter)
% interpolate missing values
% job is generally already done for seal tags (MEOP data)
% important step for float data (BGC-Argo)


%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x
% sesf01x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% last modified: 22.04.14
% 
% 20.02.11 update: addition of test to avoid error
% if no CHLA_ADJUSTED field
% if no LIGHT_ADJUSTED field
% 20.03.31 update: adding switch case for plot/no plot option (default =
% YES)
% 20.04.24 update: adjusting arrays to max_depth_dataset defined in sesf004_set_default_parameters
% 20.07.20 update: adding regularization steps for lr (low resolution)
% sealtag dataset
% 20.08.27 update: removing Chl-a values for depth > maxCHLA_depth
% 20.09.09 updates:
% reviewing fillmissing process in tag lr case
% reviewing computation of deepest measurement per profile (plot)
% 20.09.10 update: testing switch case on regulRequired_temp variable
% instead of lr/hr/fr suffix in regularization part
% 20.09.15 update: modifying displayed information at each regularixation
% step
% 20.09.23 update: transferring creation of PRES array at beginning of script
% 21.01.12 update: using CHLA instead of CHLA_ADJUSTED
% 21.03.24 update: including float data processing
% 22.03.29 update: renaming CHLA array into FLUO (Chla_nadReg becomes
% FLUO_nadReg)
% 22.04.06 update: renaming L_adxxx into PAR_nadxxx
% 22.04.14 update: TEMP and PSAL input are set to _ADJUSTED fields (contain
% density inversion correction algorithm by MEOP team)
% 
% -----------------------------------------------------------------------

disp(strcat('converting data to fixed vertical GRID for platform:',...
    platform_metadata.platform_name,'.....'))


%% set pressure array (to be used for fixed detph grid dataset)

PRES = repmat([1:max_depth_dataset]',1,np_tot) ;


%% no regularization required if platform_type = sealtag fr/hr
% (fr: high frequency, hr: "pseudo" high frequency)
% not used (20.07.20)

switch platform_metadata.platform_type
    case 'sealtag'      
        switch dataset_type(1:2)
            case 'lr'
                regulRequired_temp = 'YES' ;
            case 'hr'
                regulRequired_temp = 'NO' ;
            case 'fr'
                regulRequired_temp = 'NO' ;
        end  
    case 'float'
        regulRequired_temp = 'YES' ;
end

disp(strcat('Regularization required:',{' '},regulRequired_temp))

%% Regularisation steps if required

switch platform_metadata.platform_type
    
    
    case 'sealtag'
        
        
        switch regulRequired_temp
            
            
            case 'YES' % regularization required
                
                % pre-allocate arrays to be fit to fixed grid
                FLUO_nadReg = nan(max_depth_dataset,np_tot) ;
                PAR_nadNonLogReg = nan(max_depth_dataset,np_tot) ;
                PAR_nadLogReg = nan(max_depth_dataset,np_tot) ;
                TEMP = nan(max_depth_dataset,np_tot) ;
                SAL = nan(max_depth_dataset,np_tot) ;
                
                % interpolate FLUO/LIGHT/TEMP/PSAL
                % separate cases for TEMP and SAL / FLUO / LIGHT
                % in order to skip loop if no FLUO/LIGHT field in dataset

                
                %%%%%%%%%%%%%%%%%%%%%%%%%
                % TEMP and SAL
                %%%%%%%%%%%%%%%%%%%%%%%%%
                for ii_temp = 1:np_tot
                    disp(strcat('TEMP/SAL regularization -',...
                        ' profile nº',int2str(ii_temp),'/',int2str(np_tot)))

                    % pressure vector for profile i (irregular step)
                    pres_to_reg_temp = platform.PRES(:,ii_temp) ;

                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    % TEMP
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    % vector to be regularized
                    vect_to_reg_temp = platform.TEMP_ADJUSTED(:,ii_temp) ;
                    % output regularized vector
                    regularized_profile_temp = nan(max_depth_dataset,1) ;
                    % write in output vector where data is available
                    for jj_temp = 1:length(pres_to_reg_temp)
                        dpth_temp = pres_to_reg_temp(jj_temp) ;
                            if isnan(dpth_temp)
                            else
                                regularized_profile_temp(pres_to_reg_temp(jj_temp)) = vect_to_reg_temp(jj_temp) ;
                            end
                    end
                    % use fillmissing to overwrite nans between available data
                    regularized_profile_temp = fillmissing(regularized_profile_temp,...
                        'linear') ;
                    % first/last non nan value of profile
                    firstLastNonNanPres_temp = find(~isnan(vect_to_reg_temp)) ;
                    if isempty(firstLastNonNanPres_temp)
                        firstLastNonNanPres_temp = nan(1,2) ;
                    else
                        firstLastNonNanPres_temp = pres_to_reg_temp(firstLastNonNanPres_temp) ;
                    end
                    first_non_nan_dpth_temp = firstLastNonNanPres_temp(1) ;
                    last_non_nan_dpth_temp = firstLastNonNanPres_temp(end) ;
                    % overwrite with nans all values below last non nan depth
                    if last_non_nan_dpth_temp < max_depth_dataset
                        regularized_profile_temp(last_non_nan_dpth_temp + 1:end) = NaN ;
                    end
                    % same with all values above first non nan depth
                    % overwrite with nans all values above first non nan depth
                    if first_non_nan_dpth_temp > 1
                        regularized_profile_temp(1:first_non_nan_dpth_temp-1) = NaN ;
                    end

                    % write in data matrix
                    TEMP(:,ii_temp) = regularized_profile_temp ;
                    
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    % SAL
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    % vector to be regularized
                    vect_to_reg_temp = platform.PSAL_ADJUSTED(:,ii_temp) ;
                    % output regularized vector
                    regularized_profile_temp = nan(max_depth_dataset,1) ;
                    % write in output vector where data is available
                    for jj_temp = 1:length(pres_to_reg_temp)
                        dpth_temp = pres_to_reg_temp(jj_temp) ;
                            if isnan(dpth_temp)
                            else
                                regularized_profile_temp(pres_to_reg_temp(jj_temp)) = vect_to_reg_temp(jj_temp) ;
                            end
                    end
                    % use fillmissing to overwrite nans between available data
                    regularized_profile_temp = fillmissing(regularized_profile_temp,...
                        'linear') ;
                    % first/last non nan value of profile
                    firstLastNonNanPres_temp = find(~isnan(vect_to_reg_temp)) ;
                    if isempty(firstLastNonNanPres_temp)
                        firstLastNonNanPres_temp = nan(1,2) ;
                    else
                        firstLastNonNanPres_temp = pres_to_reg_temp(firstLastNonNanPres_temp) ;
                    end
                    first_non_nan_dpth_temp = firstLastNonNanPres_temp(1) ;
                    last_non_nan_dpth_temp = firstLastNonNanPres_temp(end) ;
                    % overwrite with nans all values below last non nan depth
                    if last_non_nan_dpth_temp < max_depth_dataset
                        regularized_profile_temp(last_non_nan_dpth_temp + 1:end) = NaN ;
                    end
                    % same with all values above first non nan depth
                    % overwrite with nans all values above first non nan depth
                    if first_non_nan_dpth_temp > 1
                        regularized_profile_temp(1:first_non_nan_dpth_temp-1) = NaN ;
                    end

                    % write in data matrix
                    SAL(:,ii_temp) = regularized_profile_temp ;
                    
                end
                
                

                %%%%%%%%%%%%%%%%%%%%%%%%%
                % CHLA_ADJUSTED
                %%%%%%%%%%%%%%%%%%%%%%%%%
                
                if isfield(platform,'CHLA_ADJUSTED') == 0
                    % no FLUO field to regularize
%                     FLUO_nadReg = nan(max_depth_dataset,np_tot) ;
                    disp('*****************************no FLUO in dataset')
                else
                    for ii_temp = 1:np_tot
                        disp(strcat('FLUO regularization -',...
                            ' profile nº',int2str(ii_temp),'/',int2str(np_tot)))

                        % pressure vector for profile i (irregular step)
                        pres_to_reg_temp = platform.PRES(:,ii_temp) ;

                        % vector to be regularized
                        vect_to_reg_temp = platform.CHLA(:,ii_temp) ;
                        % output regularized vector
                        regularized_profile_temp = nan(max_depth_dataset,1) ;
                        % write in output vector where data is available
                        for jj_temp = 1:length(pres_to_reg_temp)
                            dpth_temp = pres_to_reg_temp(jj_temp) ;
                                if isnan(dpth_temp)
                                else
                                    regularized_profile_temp(pres_to_reg_temp(jj_temp)) = vect_to_reg_temp(jj_temp) ;
                                end
                        end
                        % use fillmissing to overwrite nans between available data
                        regularized_profile_temp = fillmissing(regularized_profile_temp,...
                            'linear') ;
                        % first/last non nan value of profile
                        firstLastNonNanPres_temp = find(~isnan(vect_to_reg_temp)) ;
                        if isempty(firstLastNonNanPres_temp)
                            firstLastNonNanPres_temp = nan(1,2) ;
                        else
                            firstLastNonNanPres_temp = pres_to_reg_temp(firstLastNonNanPres_temp) ;
                        end
                        first_non_nan_dpth_temp = firstLastNonNanPres_temp(1) ;
                        last_non_nan_dpth_temp = firstLastNonNanPres_temp(end) ;
                        % overwrite with nans all values below last non nan depth
                        if last_non_nan_dpth_temp < max_depth_dataset
                            regularized_profile_temp(last_non_nan_dpth_temp + 1:end) = NaN ;
                        end
                        % same with all values above first non nan depth
                        % overwrite with nans all values above first non nan depth
                        if first_non_nan_dpth_temp > 1
                            regularized_profile_temp(1:first_non_nan_dpth_temp-1) = NaN ;
                        end

                        % write in data matrix
                        FLUO_nadReg(:,ii_temp) = regularized_profile_temp ;
                        
                    end
                    
                end
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%
                % LIGHT_ADJUSTED
                %%%%%%%%%%%%%%%%%%%%%%%%%
                
                if isfield(platform,'LIGHT_ADJUSTED') == 0
                    % no LIGHT field to regularize
%                     PAR_nadLogReg = nan(max_depth_dataset,np_tot) ;
%                     PAR_nadLogReg = nan(max_depth_dataset,np_tot) ;
                    disp('*****************************no LIGHT in dataset')
                else
                    for ii_temp = 1:np_tot
                    disp(strcat('LIGHT regularization -',...
                        ' profile nº',int2str(ii_temp),'/',int2str(np_tot)))

                        % pressure vector for profile i (irregular step)
                        pres_to_reg_temp = platform.PRES(:,ii_temp) ;

                        % vector to be regularized
                        vect_to_reg_temp = platform.LIGHT(:,ii_temp) ;
                        % output regularized vector
                        regularized_profile_temp = nan(max_depth_dataset,1) ;
                        % write in output vector where data is available
                        for jj_temp = 1:length(pres_to_reg_temp)
                            dpth_temp = pres_to_reg_temp(jj_temp) ;
                                if isnan(dpth_temp)
                                else
                                    regularized_profile_temp(pres_to_reg_temp(jj_temp)) = vect_to_reg_temp(jj_temp) ;
                                end
                        end
                        % use fillmissing to overwrite nans between available data
                        regularized_profile_temp = fillmissing(regularized_profile_temp,...
                            'linear') ;
                        % first/last non nan value of profile
                        firstLastNonNanPres_temp = find(~isnan(vect_to_reg_temp)) ;
                        if isempty(firstLastNonNanPres_temp)
                            firstLastNonNanPres_temp = nan(1,2) ;
                        else
                            firstLastNonNanPres_temp = pres_to_reg_temp(firstLastNonNanPres_temp) ;
                        end
                        first_non_nan_dpth_temp = firstLastNonNanPres_temp(1) ;
                        last_non_nan_dpth_temp = firstLastNonNanPres_temp(end) ;
                        % overwrite with nans all values below last non nan depth
                        if last_non_nan_dpth_temp < max_depth_dataset
                            regularized_profile_temp(last_non_nan_dpth_temp + 1:end) = NaN ;
                        end
                        % same with all values above first non nan depth
                        % overwrite with nans all values above first non nan depth
                        if first_non_nan_dpth_temp > 1
                            regularized_profile_temp(1:first_non_nan_dpth_temp-1) = NaN ;
                        end

                        % write in data matrix
                         PAR_nadLogReg(:,ii_temp) = regularized_profile_temp ;
                        
                    end
                    
                    % compute nonLog matrix
                    PAR_nadNonLogReg = exp(PAR_nadLogReg) ;
                    
                end


            
            otherwise % platform dataset should be fr/hr
                

                disp('platform_type: sealtag fr/hr => no regularization required')
                % FLUO
                if isfield(platform,'CHLA_ADJUSTED')
                    FLUO_nadReg = platform.CHLA(1:max_depth_dataset,:) ;
                else
                    FLUO_nadReg = nan(max_depth_dataset,np_tot) ;
                    disp('*****************************no FLUO in dataset')
                end
                % LIGHT
                if isfield(platform,'LIGHT_ADJUSTED')
                    PAR_nadLogReg = platform.LIGHT_ADJUSTED(1:max_depth_dataset,:) ;
                    PAR_nadNonLogReg = exp(platform.LIGHT_ADJUSTED(1:max_depth_dataset,:)) ;
                else
                    PAR_nadLogReg = nan(max_depth_dataset,np_tot) ;
                    PAR_nadNonLogReg = nan(max_depth_dataset,np_tot) ;
                    disp('*****************************no LIGHT in dataset')
                end
                
                % TEMP & SAL
                TEMP = platform.TEMP_ADJUSTED(1:max_depth_dataset,:) ;
                SAL = platform.PSAL_ADJUSTED(1:max_depth_dataset,:) ;
                
        end

    
    case 'float'
        %% set arrays to be fit to fixed grid

        FLUO_nadReg = nan(max_depth_dataset,np_tot) ;
        PAR_nadNonLogReg = nan(max_depth_dataset,np_tot) ;
        PAR_nadLogReg = nan(max_depth_dataset,np_tot) ;
        TEMP = nan(max_depth_dataset,np_tot) ;
        SAL = nan(max_depth_dataset,np_tot) ;




        %% interpolate TEMP/PSAL/FLUO/LIGHT

        for ii_temp = 1:np_tot
            disp(strcat('TEMP/SAL/FLUO/LIGHT regularization -',...
                ' profile nº',int2str(ii_temp),'/',int2str(np_tot)))

            % pressure vector for profile i (irregular step)
            pres_to_reg_temp = platform.PRES(:,ii_temp) ;
            % leave all parameters set at nan if no measurement below 2
            % m depth (i.e. at least 2 measurements in regularized profile,
            % at 1 m and at 2 m)
            if max(pres_to_reg_temp) < 2
            else

                %%%%%%%%%%%%%%%%%%%%%%%%%
                % TEMP
                %%%%%%%%%%%%%%%%%%%%%%%%%
                % vector to be regularized
                vect_to_reg_temp = platform.TEMP_ADJUSTED(:,ii_temp) ;
                % output regularized vector
                regularized_profile_temp = nan(max_depth_dataset,1) ;

                % number of non nan values in vect_to_reg_temp
                n_values_temp = nnz(~isnan(vect_to_reg_temp)) ;

                % leave all nan in xxxReg vector if only 1 non nan value in raw dataset
                if n_values_temp < 2
                else

                    % replace nans in vector to be interpolated
                    newVect_to_reg_temp = fillmissing(vect_to_reg_temp,'linear') ;

                    % index of last non NaN value
                    bot_temp = min(find(~isnan(vect_to_reg_temp),1,'last')) ;
                    % index of first non NaN value
                    top_temp = max(find(~isnan(vect_to_reg_temp),1,'first')) ;
                    % first non NaN and last non NaN in meters
                    botm_temp = floor(pres_to_reg_temp(bot_temp)) ;
                    topm_temp = max(1,ceil(pres_to_reg_temp(top_temp))) ;

                    % consolidate (remove duplicates)
                    % consolidator function by by John D'Errico
                    % https://www.mathworks.com/matlabcentral/fileexchange/8354-consolidator
                    X_temp = pres_to_reg_temp(top_temp:bot_temp) ;
                    Y_temp = newVect_to_reg_temp(top_temp:bot_temp) ;
                    [xg_temp,yg_temp] = consolidator(X_temp,Y_temp,'mean') ;

                    % write interpolated values according to 1m delta PRES array
                    interp_vect_temp = interp1(xg_temp,yg_temp,[topm_temp:botm_temp]) ;

                    % write interpolated values in xxxReg(:,i) vector of profile i
                    % within interpolated array bounds
                    regularized_profile_temp(topm_temp:botm_temp) =...
                        interp_vect_temp ;
                    TEMP(:,ii_temp) = regularized_profile_temp ;
                end


                %%%%%%%%%%%%%%%%%%%%%%%%%
                % SAL
                %%%%%%%%%%%%%%%%%%%%%%%%%
                % vector to be regularized
                vect_to_reg_temp = platform.PSAL_ADJUSTED(:,ii_temp) ;
                % output regularized vector
                regularized_profile_temp = nan(max_depth_dataset,1) ;

                % number of non nan values in vect_to_reg_temp
                n_values_temp = nnz(~isnan(vect_to_reg_temp)) ;

                % leave all nan in xxxReg vector if only 1 non nan value in raw dataset
                if n_values_temp < 2
                else

                    % replace nans in vector to be interpolated
                    newVect_to_reg_temp = fillmissing(vect_to_reg_temp,'linear') ;

                    % index of last non NaN value
                    bot_temp = min(...
                        find(~isnan(vect_to_reg_temp),1,'last'),...
                        find(pres_to_reg_temp <= max_depth_dataset,1,'last')) ;
                    % index of first non NaN value
                    top_temp = max(...
                        find(~isnan(vect_to_reg_temp),1,'first'),...
                        find(pres_to_reg_temp > 1,1,'first')) ;
                    % (additinal conditions are set for bot_temp top/temp
                    % because max_depth_dataset is computed with parameter TEMP)

                    % first non NaN and last non NaN in meters
                    botm_temp = floor(pres_to_reg_temp(bot_temp)) ;
                    topm_temp = floor(pres_to_reg_temp(top_temp)) ;

                    % consolidate (remove duplicates)
                    % consolidator function by by John D'Errico
                    % https://www.mathworks.com/matlabcentral/fileexchange/8354-consolidator
                    X_temp = pres_to_reg_temp(top_temp:bot_temp) ;
                    Y_temp = newVect_to_reg_temp(top_temp:bot_temp) ;
                    [xg_temp,yg_temp] = consolidator(X_temp,Y_temp,'mean') ;

                    % write interpolated values according to 1m delta PRES array
                    interp_vect_temp = interp1(xg_temp,yg_temp,[topm_temp:botm_temp]) ;

                    % write interpolated values in xxxReg(:,i) vector of profile i
                    % within interpolated array bounds
                    regularized_profile_temp(topm_temp:botm_temp) =...
                        interp_vect_temp ;
                    SAL(:,ii_temp) = regularized_profile_temp ;
                end


                %%%%%%%%%%%%%%%%%%%%%%%%%
                % FLUO
                %%%%%%%%%%%%%%%%%%%%%%%%%
                % vector to be regularized
                vect_to_reg_temp = platform.CHLA(:,ii_temp) ;
                % output regularized vector
                regularized_profile_temp = nan(max_depth_dataset,1) ;

                % number of non nan values in vect_to_reg_temp
                n_values_temp = nnz(~isnan(vect_to_reg_temp)) ;

                % leave all nan in xxxReg vector if only 1 non nan value in raw dataset
                if n_values_temp < 2
                else

                    % replace nans in vector to be interpolated
                    newVect_to_reg_temp = fillmissing(vect_to_reg_temp,'linear') ;

                    % index of last non NaN value
                    bot_temp = min(...
                        find(~isnan(vect_to_reg_temp),1,'last'),...
                        find(pres_to_reg_temp <= max_depth_dataset,1,'last')) ;
                    % index of first non NaN value
                    top_temp = max(...
                        find(~isnan(vect_to_reg_temp),1,'first'),...
                        find(pres_to_reg_temp > 1,1,'first')) ;
                    % (additinal conditions are set for bot_temp top/temp
                    % because max_depth_dataset is computed with parameter TEMP)
                    % first non NaN and last non NaN in meters
                    botm_temp = floor(pres_to_reg_temp(bot_temp)) ;
                    topm_temp = floor(pres_to_reg_temp(top_temp)) ;

                    % consolidate (remove duplicates)
                    % consolidator function by by John D'Errico
                    % https://www.mathworks.com/matlabcentral/fileexchange/8354-consolidator
                    X_temp = pres_to_reg_temp(top_temp:bot_temp) ;
                    Y_temp = newVect_to_reg_temp(top_temp:bot_temp) ;
                    [xg_temp,yg_temp] = consolidator(X_temp,Y_temp,'mean') ;

                    % write interpolated values according to 1m delta PRES array
                    interp_vect_temp = interp1(xg_temp,yg_temp,[topm_temp:botm_temp]) ;

                    % write interpolated values in xxxReg(:,i) vector of profile i
                    % within interpolated array bounds
                    regularized_profile_temp(topm_temp:botm_temp) =...
                        interp_vect_temp ;
                    FLUO_nadReg(:,ii_temp) = regularized_profile_temp ;
                end


                %%%%%%%%%%%%%%%%%%%%%%%%%
                % LIGHT
                %%%%%%%%%%%%%%%%%%%%%%%%%
                % vector to be regularized
                vect_to_reg_temp = platform.DOWNWELLING_PAR(:,ii_temp) ;
                % output regularized vector
                regularized_profile_temp = nan(max_depth_dataset,1) ;

                % number of non nan values in vect_to_reg_temp
                n_values_temp = nnz(~isnan(vect_to_reg_temp)) ;

                % leave all nan in xxxReg vector if only 1 non nan value in raw dataset
                if n_values_temp < 2
                else

                    % replace nans in vector to be interpolated
                    newVect_to_reg_temp = fillmissing(vect_to_reg_temp,'linear') ;

                    % index of last non NaN value
                    bot_temp = min(...
                        find(~isnan(vect_to_reg_temp),1,'last'),...
                        find(pres_to_reg_temp <= max_depth_dataset,1,'last')) ;
                    % index of first non NaN value
                    top_temp = max(...
                        find(~isnan(vect_to_reg_temp),1,'first'),...
                        find(pres_to_reg_temp > 1,1,'first')) ;
                    % (additinal conditions are set for bot_temp top/temp
                    % because max_depth_dataset is computed with parameter TEMP)
                    % first non NaN and last non NaN in meters
                    botm_temp = floor(pres_to_reg_temp(bot_temp)) ;
                    topm_temp = floor(pres_to_reg_temp(top_temp)) ;

                    % consolidate (remove duplicates)
                    % consolidator function by by John D'Errico
                    % https://www.mathworks.com/matlabcentral/fileexchange/8354-consolidator
                    X_temp = pres_to_reg_temp(top_temp:bot_temp) ;
                    Y_temp = newVect_to_reg_temp(top_temp:bot_temp) ;
                    [xg_temp,yg_temp] = consolidator(X_temp,Y_temp,'mean') ;

                    % write interpolated values according to 1m delta PRES array
                    interp_vect_temp = interp1(xg_temp,yg_temp,[topm_temp:botm_temp]) ;

                    % write interpolated values in xxxReg(:,i) vector of profile i
                    % within interpolated array bounds
                    regularized_profile_temp(topm_temp:botm_temp) =...
                        interp_vect_temp ;
                    PAR_nadNonLogReg(:,ii_temp) = regularized_profile_temp ;
                end
                
            end

        end
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional plot
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % plot raw and reg values
            %%%%%%%%%%%%%%%%%%%%%%%%% 
        %     figure(6021)
        %     sgtitle(strcat('profile',int2str(i)))
        %     
        %     % FLUO
        %     subplot(1,4,1)
        %     plot(FLUO_nadReg(:,i),PRES(:,i),'-g')
        %     hold on  
        %     plot(platform.CHLA(:,i),platform.PRES(:,i),'xk')
        %     plot(FLUO_nadReg(:,i),PRES(:,i),'.k')
        %     hold off
        %     ylim([0 200])
        %     legend({'reg profile','raw data','added points'})
        %     
        % 	% LIGHT
        %     subplot(1,4,2)
        %     plot(PAR_nadLogReg(:,i),PRES(:,i),'-g')
        %     hold on  
        %     plot(platform.DOWNWELLING_PAR(:,i),platform.PRES(:,i),'xk')
        %     plot(PAR_nadLogReg(:,i),PRES(:,i),'.k')
        %     hold off
        %     ylim([0 200])
        %     legend({'reg profile','raw data','added points'})
        %     
        %     % TEMP
        %     subplot(1,4,3)
        %     plot(TEMP(:,i),PRES(:,i),'-g')
        %     hold on  
        %     plot(platform.TEMP_ADJUSTED(:,i),platform.PRES(:,i),'xk')
        %     plot(TEMP(:,i),PRES(:,i),'.k')
        %     hold off
        %     ylim([0 200])
        %     legend({'reg profile','raw data','added points'})
        %     
        %     % SAL
        %     subplot(1,4,4)
        %     plot(SAL(:,i),PRES(:,i),'-g')
        %     hold on  
        %     plot(platform.PSAL_ADJUSTED(:,i),platform.PRES(:,i),'xk')
        %     plot(SAL(:,i),PRES(:,i),'.k')
        %     hold off
        %     ylim([0 200])
        %     legend({'reg profile','raw data','added points'})



end


%% plot time series of max depth per profile

switch doPlots
    case 'NO'
        
    case 'YES'
        deepestMeasPerProfile_temp = find_ndim(~isnan(TEMP),1,'last') ;
        [maxdepth_temp,indmaxdepth_temp] = max(deepestMeasPerProfile_temp) ;
        figure (121)
        plot1_temp = plot([1:np_tot],deepestMeasPerProfile_temp,'.k',...
            'DisplayName','deepest measurement per profile') ;
        hold on
        plot2_temp = plot(indmaxdepth_temp, maxdepth_temp, 'or',...
            'DisplayName',...
            strcat('deepest measurement in dataset=',int2str(maxdepth_temp),' m')) ;
        title(strcat('max depth per profile - platform:',platform_metadata.platform_name))
        set(gca,'YDir','reverse')
        ylim([-100 max_depth_dataset])
        xlabel('profile nº')
        ylabel('depth (m)')
        legend([plot2_temp], 'Location','NorthEast')
        hold off

end



%% remove values/insert NaN according to user-defined max_depth variable

diff_temp = max_depth - max_depth_dataset ;
if diff_temp > 0
    disp(strcat('user-defined max_depth (', int2str(max_depth),...
        ' m) is deeper than dataset maximum depth (', int2str(max_depth_dataset),...
        ' m)'))
    disp('====> bottom part of dataset profiles is filled with NaN values')
    intervalToFill_temp = max_depth_dataset + 1:max_depth ;
    PRES(intervalToFill_temp,:) = NaN ;
    FLUO_nadReg(intervalToFill_temp,:) = NaN ;
    PAR_nadNonLogReg(intervalToFill_temp,:) = NaN ;
    PAR_nadLogReg(intervalToFill_temp,:) = NaN ;
    TEMP(intervalToFill_temp,:) = NaN ;
    SAL(intervalToFill_temp,:) = NaN ;
elseif diff_temp < 0
    disp(strcat('user-defined max_depth (', int2str(max_depth),...
        ' m) is shallower than dataset maximum depth (', int2str(max_depth_dataset),...
        ' m)'))
    disp(char(strcat('====> bottom part of profiles deeper than',{' '},int2str(max_depth),...
        ' m are cut')))
    intervalToFill_temp = max_depth + 1:max_depth_dataset ;
    PRES(intervalToFill_temp,:) = [] ;
    FLUO_nadReg(intervalToFill_temp,:) = [] ;
    PAR_nadNonLogReg(intervalToFill_temp,:) = [] ;
    PAR_nadLogReg(intervalToFill_temp,:) = [] ;
    TEMP(intervalToFill_temp,:) = [] ;
    SAL(intervalToFill_temp,:) = [] ;
elseif diff_temp == 0
end

    


%% clear temp variables
clear *_temp

%% END