%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% preprocess PAR data
% compute Zeu euphotic depth / dquench quenching depth
% Zpd/Zeu/Quenching values are computed on PAR profile functional fit (strictly decreasing with depth)

%%%%%%%%%% REFERENCES %%%%%%%%%%
% Zeu: Morel and Berthon, 1989
% Zpd:
% Howard & McCluney, 1975 https://doi.org/10.1364/AO.14.000413
% Morel, 1988 https://doi.org/10.1029/JC093iC09p10749
% quenching: Xing et al. 2018

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
% J. Weis (jw)
% last modified: 22.05.10
% 
% 21.04.27 update: gathering metrics in same script
% 22.02.03 updates:
% including saturationDepth in subsurValues computation
% adding ZeuInterp variable to lightData table (interp. window: approx. 1
% day)
% computing Zpd from ZeuInterp instead of Zeu "raw"
% 22.04.06 update: renaming LIGHT into PAR (L_adxxx becomes PAR_nadxxx)
% 22.04.18 update: adding computation of mean Kd
% 22.05.10 update: computation of meanKd is made on functional fit + adding
% computation of meanKdML the mean Kd in the mixed layer
% 
% -----------------------------------------------------------------------

disp(strcat('PAR data pre-processing step7: post-processing METRICS / platform:',...
    platform_metadata.platform_name,'.....'))


%% new subsurface value (functional fit)
% subsurIndexes_temp = sub2ind(size(PAR_nadLogRegDkSaFitAll),...
%     parData.saturationDepth,...
%     transpose(1:platform_metadata.np)) ;
% subsurValues_temp = nan(platform_metadata.np,1) ;
% subsurValues_temp(~isnan(subsurIndexes_temp)) =...
%     PAR_nadLogRegDkSaFitAll(subsurIndexes_temp(~isnan(subsurIndexes_temp))) ;
% parData.subsurValFit = subsurValues_temp ;
% +++++++++++ FIXING SUB2IND ISSUE ++++++++++++++
finiteIndices = find(isfinite(parData.saturationDepth));
parData.subsurValFit = NaN(platform_metadata.np,1);
parData.subsurValFit(finiteIndices,1) = arrayfun(@(x) PAR_nadLogRegDkSaFitAll(parData.saturationDepth(x), x), finiteIndices);


%% compute euphotic depth
% defined as depth where the downwelling PAR irradiance
% is reduced to 1% of its value at the surface
% Morel, A., and J.-F. Berthon (1989)

zeuDep_temp = nan(platform_metadata.np,1) ;
zeuzeuDep_temp = arrayfun(@(a) find(PAR_nadNonLogRegDkSaFitAll(:,a) / exp(parData.subsurValFit(a)) < 0.01,...
    1,'first'),find(idx031_lightNonNan),...
    'UniformOutput',false) ;
zeuNonEmptyIndex_temp = arrayfun(@(a) ~isempty(zeuzeuDep_temp{a,1}),...
    [1:numel(find(idx031_lightNonNan))]) ;
zeuDep_temp(zeuNonEmptyIndex_temp) = cell2mat(zeuzeuDep_temp(zeuNonEmptyIndex_temp)) ;

parData.Zeu = zeuDep_temp ;

% fill missing and smooth Zeu data (window size: approx. 1 day)
zeuDepFilled_temp  = fillmissing(zeuDep_temp,...
    'movmedian',platform_metadata.np / max(genData.depDayNo),...
    'EndValues','none') ;
parData.ZeuInterp = movmedian(zeuDepFilled_temp,platform_metadata.np / max(genData.depDayNo),...
    'Endpoints','shrink') ;


%% penetration depth (first optical depth)

parData.Zpd = parData.ZeuInterp / 4.6 ;

%% compute depth of quenching threshold 15 umol.m-2.s-1 (xing et al. 2018)

% parData.quenchDepth contains depth of quenching threshold at 15 umol.m-2.s-1
% (xing et al. 2018)

quenchDep_temp = find_ndim(PAR_nadNonLogRegDkSaFitAll > 15, 1, 'last').' ;
quenchDep_temp(quenchDep_temp == 0) = NaN ;

parData.quenchDepth = quenchDep_temp ;

%% compute mean Kd / mean Kd in Mixed Layer (ML)

% mean Kd
meanKd_temp = arrayfun(@(a) mean(...
    - diff(PAR_nadLogRegDkSaFitAll(max(1,parData.saturationDepth(a)):min(parData.lastNonNan(a),1000),a)),...
    'omitnan'),...
    1:length(platform.LATITUDE)).' ;

parData.meanKd = meanKd_temp ;

% mean Kd ML
meanKdML_temp = arrayfun(@(a) mean(...
    - diff(PAR_nadLogRegDkSaFitAll(max(1,parData.saturationDepth(a)):min(genData.MLDphy(a),1000),a)),...
    'omitnan'),...
    1:length(platform.LATITUDE)).' ;

parData.meanKdML = meanKdML_temp ;



%% clear temp data
clear *_temp

%% END