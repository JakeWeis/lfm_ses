%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% compute profile metrics after basic data processing
% (dark + NPQ + smoothing/functional fit)
% metrics include:
% surfFluo / maxFluo / depth of maxFluo / intFluo / MLDbio /
% smoothing residuals (%)

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x
% sesf01x
% sesf03x
% sesf04x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% last modified: 22.03.29
% 
% -----------------------------------------------------------------------

disp(strcat('FLUO data pre-processing step7: post-processing METRICS / platform:',...
    platform_metadata.platform_name,'.....'))


%% surf
% use of FLUO_nadRegDkNpq

surfFluoVal_temp = FLUO_nadRegDkNpq(1,:).' ;

%% surf as mean in first optical depth (penetration depth Zpd) - for satellite matchup

zpd_temp = floor(fillmissing(parData.Zpd,'nearest')) ;
meanFluoZpd_temp = nan(platform_metadata.np,1) ;
meanFluoZpd_temp(idx041_chlaNonNan) =...
    arrayfun(@(a) mean(FLUO_nadRegDkNpq(1:zpd_temp(a),a),'omitnan'),...
    find(idx041_chlaNonNan)) ;

%% mean Fluo in MLD
zmld_temp = fillmissing(genData.MLDphy,'nearest') ;
meanFluoZmld_temp = nan(platform_metadata.np,1) ;
meanFluoZmld_temp(idx041_chlaNonNan) =...
    arrayfun(@(a) mean(FLUO_nadRegDkNpqFitBnd(1:zmld_temp(a),a),'omitnan'),...
    find(idx041_chlaNonNan)) ;


%% max/depth of max -> based on smoothed profile
% use of FLUO_nadRegDkNpq + FLUO_nadRegDkNpqSth

% depth of max -> based on smoothed profile
[maxOfFFit_temp,maxFluoDepth_temp] = max(FLUO_nadRegDkNpqFitAll,[],1) ;
% unsmoothed max value of Chl-a -> retrieve corresponding unsmoothed value
maxFluoIndexes_temp = sub2ind(size(FLUO_nadRegDkNpq),...
    maxFluoDepth_temp.',...
    transpose(1:platform_metadata.np)) ;
maxFluoValues_temp = nan(platform_metadata.np,1) ;
maxFluoValues_temp(~isnan(maxFluoIndexes_temp)) =...
    FLUO_nadRegDkNpq(maxFluoIndexes_temp(~isnan(maxFluoIndexes_temp))) ;
maxFluoDepth_temp(~idx041_chlaNonNan) = NaN ;



%% int FLUO
% use of FLUO_nadRegDkNpqFitAll

% initialize integFluo_temp
integFLuo_temp = nan(np_tot,1) ;
% use modified FLUO_nadRegDkNpq array so as not to have NaN values
% hypothesis: zero values do have very few contribution to integration)
fluoForInteg_temp = FLUO_nadRegDkNpqFitAll ;
fluoForInteg_temp(isnan(fluoForInteg_temp)) = 0 ;
% vertical vector of depths to be used in integration
presForInteg_temp = [1:max_depth] ;
% compute trapezoidal integration
integFLuo_temp = trapz(presForInteg_temp,fluoForInteg_temp,1) ;


%% MLDbio
%%%%% METHOD1 %%%%%
% % use of FLUO_nadRegDkNpqSth
% 
% % initialize variables
% dfluo_dpres = nan(max_depth - 1,np_tot) ;
% % vertical length of dfluo_dpres is (max_depth - 1)
% % because of the further use of diff function
% MLDbio_temp = nan(np_tot,1) ;
% 
% for ii_temp = find(idx041_chlaNonNan).'
%     % compute smoothed fluo local derivative
%     dfluo_dpres(:,ii_temp) =...
%         diff(FLUO_nadRegDkNpqSth(:,ii_temp))./diff(PRES(:,ii_temp)) ;
%     % clear discontinuities of the derivative
%     % other method: apply additonal smoothing on FLUO_nadRegDkNpqSth
%     indinf_temp = find(dfluo_dpres(:,ii_temp) == -Inf) ;
%     dfluo_dpres(indinf_temp,ii_temp) = 0 ;
%     % get index and depth of max slope (negative slope => min dfluo/dpres)
%     % use of shift_temp for data consistency: smoothing filter has a
%     % 5-point + 7-point window
%     % => derivative makes sense for first non NaN depth + 5 meters
%     shift_temp = fluoData.firstNonNan(ii_temp) + 5 ; % adding 5m margin
%     [M_temp,mldb_temp] = min(dfluo_dpres(shift_temp:end,ii_temp)) ;
%     % correct the computed index with index shift
%     mldb_temp = mldb_temp + shift_temp - 1 ;
%     % write in MLDbio array
%     MLDbio_temp(ii_temp) = PRES(mldb_temp,ii_temp) ;
% end


%%%%% METHOD2 %%%%%
% use of FLUO_nadRegDkNpqFitAll
indexes_temp =...
    idx000_genProfileSelec &...
    idx041_chlaNonNan ;  
indFitAll_temp = find(indexes_temp) ;

MLDbio_temp = nan(np_tot,1) ;

for iProfile_temp = indFitAll_temp.'
    coef_temp = fluo_fdAllCoefs(:,iProfile_temp) ;
    vec0_temp = FLUO_nadRegDkNpqFitAll(:,iProfile_temp) ;
    pp_temp = find(~isnan(vec0_temp)) ;
    wbasis_temp = create_bspline_basis([pp_temp(1) pp_temp(end)],nbaz,nord) ;
    objfd_temp = fd(coef_temp,wbasis_temp) ;
    dobjdpresfd_temp = deriv_fd(objfd_temp) ;
    dobjdpres_temp = eval_fd(dobjdpresfd_temp,pp_temp) ;
    [dmin_temp,indMin_temp] = min(dobjdpres_temp) ;
    MLDbio_temp(iProfile_temp) = pp_temp(indMin_temp) ;
end


%% smoothing residuals relative to smoothed profile
% use of FLUO_nadRegDkNpq + use of FLUO_nadRegDkNpqSth

% initialize smoothRelResid_temp
smoothRelResid_temp = nan(np_tot,1) ;

smoothResid_temp = abs(FLUO_nadRegDkNpq-FLUO_nadRegDkNpqFitAll) ;
smoothRelResid_temp = smoothResid_temp ./ FLUO_nadRegDkNpqFitAll ;
smoothRelResid_temp = median(smoothRelResid_temp,1,'omitnan') ;


%% write data in fluoData table
fluoData.surf = surfFluoVal_temp ;
fluoData.surfZpd = meanFluoZpd_temp ;
fluoData.meanFluoMLD = meanFluoZmld_temp ;
fluoData.maxVal = maxFluoValues_temp ;
fluoData.maxValFit = maxOfFFit_temp.' ;
fluoData.maxValDepth = maxFluoDepth_temp.' ;
fluoData.integ = integFLuo_temp.' ;
fluoData.MLDbio = MLDbio_temp ;
fluoData.smoothRelResid = smoothRelResid_temp.' ;

% write nan where no Chl-a profiles available
fluoData(~idx041_chlaNonNan,2:end) = array2table(...
    nan(nnz(~idx041_chlaNonNan),...
numel(fluoData.Properties.VariableNames) - 1)) ;


%% clear temp variables
clear *_temp

%% END
