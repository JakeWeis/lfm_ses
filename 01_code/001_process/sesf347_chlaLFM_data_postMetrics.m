%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% compute ChlaLFM profile metrics after LFM prediction
% metrics include:
% surfChla / maxChla / depth of maxChla / intChla

%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x
% sesf01x
% sesf03x
% sesf04x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
disp('***MATLAB version information***')
disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% last modified: 22.10.28
% 
% 22.03.31 update: add creation of chlaData table (based on fluoData table)
% 22.10.28 update: converting script for CHLA_LFM
% 
% -----------------------------------------------------------------------

disp(strcat('CHLA_LFM data post-processing METRICS / platform:',...
    platform_metadata.platform_name,'.....'))


%% surf
% use of CHLA_adRegDkNpqOc

surfChlaVal_temp = CHLA_LFM(1,:).' ;

%% surf as mean in first optical depth (penetration depth Zpd) - for satellite matchup

zpd_temp = floor(fillmissing(parData.Zpd,'nearest')) ;
meanChlaZpd_temp = nan(platform_metadata.np,1) ;
meanChlaZpd_temp(idx036_lightFDAfitBnd) =...
    arrayfun(@(a) mean(CHLA_LFM(1:zpd_temp(a),a),'omitnan'),...
    find(idx036_lightFDAfitBnd)) ;

%% mean Chla in MLD
zmld_temp = fillmissing(genData.MLDphy,'nearest') ;
meanChlaZmld_temp = nan(platform_metadata.np,1) ;
meanChlaZmld_temp(idx036_lightFDAfitBnd) =...
    arrayfun(@(a) mean(CHLA_LFM(1:zmld_temp(a),a),'omitnan'),...
    find(idx036_lightFDAfitBnd)) ;


%% max/depth of max -> based on smoothed profile
% use of CHLA_adRegDkNpqOc + CHLA_adRegDkNpqOcFitAll

% depth of max -> based on smoothed profile
[maxOfFFit_temp,maxChlaDepth_temp] = max(CHLA_LFM,[],1) ;
% unsmoothed max value of Chl-a -> retrieve corresponding unsmoothed value
maxChlaIndexes_temp = sub2ind(size(CHLA_LFM),...
    maxChlaDepth_temp.',...
    transpose(1:platform_metadata.np)) ;
maxChlaValues_temp = nan(platform_metadata.np,1) ;
maxChlaValues_temp(~isnan(maxChlaIndexes_temp)) =...
    CHLA_LFM(maxChlaIndexes_temp(~isnan(maxChlaIndexes_temp))) ;
maxChlaDepth_temp(~idx036_lightFDAfitBnd) = NaN ;



%% int Chla
% use of CHLA_adRegDkNpqOcFitAll

% initialize integChla_temp
integChla_temp = nan(platform_metadata.np,1) ;
% use modified CHLA array so as not to have NaN values - see hereafter
% hypothesis: zero values do have very few contribution to integration)
chlaForInteg_temp = CHLA_LFM ;
chlaForInteg_temp(isnan(chlaForInteg_temp)) = 0 ;
% vertical vector of depths to be used in integration
presForInteg_temp = [1:max_depth] ;
% compute trapezoidal integration
integChla_temp = trapz(presForInteg_temp,chlaForInteg_temp,1) ;


%% MLDbio

MLDbio_temp = nan(platform_metadata.np,1) ;


%% smoothing residuals relative to smoothed profile
% use of CHLA_adRegDkNpqOc + use of CHLA_adRegDkNpqOcSth

% initialize smoothRelResid_temp
smoothRelResid_temp = nan(platform_metadata.np,1) ;


%% write data in chlaData table

% create chlaData table
chlaLFMData = table() ;
chlaLFMData.profileNB = fluoData.profileNB ;

% write data
chlaLFMData.surf = surfChlaVal_temp ;
chlaLFMData.surfZpd = meanChlaZpd_temp ;
chlaLFMData.meanChlaMLD = meanChlaZmld_temp ;
chlaLFMData.maxVal = maxChlaValues_temp ;
chlaLFMData.maxValFit = maxOfFFit_temp.' ;
chlaLFMData.maxValDepth = maxChlaDepth_temp.' ;
chlaLFMData.integ = integChla_temp.' ;
chlaLFMData.MLDbio = MLDbio_temp ;
chlaLFMData.smoothRelResid = smoothRelResid_temp ;

% write nan where no Chl-a profiles available
chlaLFMData(~idx036_lightFDAfitBnd,2:end) = array2table(...
    nan(nnz(~idx036_lightFDAfitBnd),...
numel(chlaLFMData.Properties.VariableNames) - 1)) ;


%% clear temp variables
clear *_temp

%% END
