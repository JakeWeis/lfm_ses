%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% preprocess PAR data
% proceed to functional fit
% 2 fits are computed:
%     - 1 fit on a fixed interval (min_topIntegBound:min_botIntegBound)
%     - 1 fit on the entire vertical profile where data is non Nan
% Compute non-log PAR from functional fit (computed on log(PAR))

%%%%%%%%%% REFERENCES %%%%%%%%%%
% fit with constraints from fda package
% https://www.psych.mcgill.ca/misc/fda/software.html

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
% disp('***MATLAB version information***')
% disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
% disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% J. Weis (jw)
% last modified: 22.04.06
% 
% 22.02.04 update: adding idxARGVALS condition on profile length
% (ARGVALS in smooth_monotone needs at least three values)
% 22.04.06 update: renaming LIGHT into PAR (L_adxxx becomes PAR_nadxxx)
% 22.10.28 update: adding computation of Kd from PAR
% 
% ----------------------------------------------------------------------

disp(strcat('PAR data pre-processing step6: fda fit / platform:',...
    platform_metadata.platform_name,'.....'))


%% INITIALIZE DATA

% fitted data (entire interval)
parFdaFitAll_temp = nan(max_depth,platform_metadata.np) ;
PAR_nadLogRegDkSaFitAll = nan(max_depth,platform_metadata.np) ;

% fitted data (fixed interval)
parFdaFitBnd_temp = nan(max_depth,platform_metadata.np) ;
PAR_nadLogRegDkSaFitBnd = nan(max_depth,platform_metadata.np) ;
Kd_Bnd = nan(max_depth,platform_metadata.np) ;

% fd objects
par_fdAllCoefs = nan(nbaz,platform_metadata.np) ;
par_fdBndCoefs = nan(nbaz,platform_metadata.np) ;
KdFDcoefs = nan(nbaz,platform_metadata.np) ;

%% COMPUTE FIT ON ENTIRE VERTICAL PROFILE
% SELECT PROFILES TO COMPUTE
lumToFitAll_temp = PAR_nadLogRegDkSa ;

% ARGVALS in smooth_monotone needs at least three values.
idxARGVALS_temp = sum(~isnan(lumToFitAll_temp)).' >= 3 ;
indexes_temp =...
    idx000_genProfileSelec &...
    idx031_lightNonNan &...
    idx035_lightDay &...
    idxARGVALS_temp ;
indFitAll_temp = find(indexes_temp) ;

% DATA2FD PARAMETERS
% define basis
nbaz_temp = nbaz ;
nord_temp = nord ;

% FIT PAR
%  starting values for coefficient
cvecAll_temp = zeros(nbaz_temp,1) ;
%  set up functional parameter object
Lfdobj_temp    = 1  ;    %  original value in script = 2 ; penalize curvature of acceleration
lambda_temp    = lambda  ; %  smoothing parameter %10^(-0.1)
% compute monotone fit
for iProfile_temp = indFitAll_temp.'
    % display computation progress every 100
    if rem(iProfile_temp,100) == 0
        disp(strcat('ffit platform:',32,char(platform_metadata.platform_name),...
            '/profile:',32,int2str(iProfile_temp),'/',int2str(platform_metadata.np)))
    end
    % parameters
    vec0_temp = lumToFitAll_temp(:,iProfile_temp) ;
    pp_temp = find(~isnan(vec0_temp)) ;
    if ~isempty(pp_temp)
        vec_temp = vec0_temp(pp_temp) ;
        wbasis_temp = create_bspline_basis([pp_temp(1) pp_temp(end)],nbaz_temp,nord_temp) ;
        Wfd0_temp  = fd(cvecAll_temp, wbasis_temp) ;
        proffdPar_temp = fdPar(Wfd0_temp,Lfdobj_temp,lambda_temp) ;
        % monotone smooth
        [~, ~, vecfd0_temp, ~, ~, ~, ~] =...
            smooth_monotone(pp_temp,vec_temp,proffdPar_temp) ;
        % eval fit
        evalvecfd_temp = eval_fd(vecfd0_temp,pp_temp) ;
        parFdaFitAll_temp(pp_temp,iProfile_temp) = evalvecfd_temp ;
        par_fdAllCoefs(:,iProfile_temp) = getcoef(vecfd0_temp) ;
    end
end

% write data
PAR_nadLogRegDkSaFitAll = parFdaFitAll_temp ;


%% COMPUTE FIT ON A FIXED INTERVAL
% SELECT PROFILES TO COMPUTE
lumToFitBnd_temp = PAR_nadLogRegDkSa ;
idxBnd_temp =...
    find_ndim(~isnan(lumToFitBnd_temp),1,'first') <= min_topIntegBound &...
    find_ndim(~isnan(lumToFitBnd_temp),1,'last') >= min_botIntegBound ;
idx036_lightFDAfitBnd =...
    idxBnd_temp.' &...
    idx000_genProfileSelec &...
    idx035_lightDay ;
indFitBnd_temp = find(idx036_lightFDAfitBnd) ;

% DATA2FD PARAMETERS
% depth interval
pp_temp = (min_topIntegBound:1:min_botIntegBound).' ;
% Data to fit/convert to fd object
lumToFdBnd_temp = lumToFitBnd_temp(pp_temp,indFitBnd_temp) ;
% define basis
nbaz_temp = nbaz ;
nord_temp = nord ;
wbasis_temp = mybL ;
% set up functional parameter object
cvecBnd_temp = zeros(nbaz_temp,length(indFitBnd_temp)) ;
Lfdobj_temp    = 1  ;    %  original value in script = 2 ; penalize curvature of acceleration
lambda_temp    = 0.08  ; %  smoothing parameter %10^(-0.1)
Wfd0_temp  = fd(cvecBnd_temp, wbasis_temp) ;
proffdPar_temp = fdPar(Wfd0_temp,Lfdobj_temp,lambda_temp) ;

% compute monotone fit
[~, ~, vecfdBnd_temp, ~, ~, ~, ~] =...
    smooth_monotone(pp_temp,lumToFdBnd_temp,proffdPar_temp) ;
% EVAL FIT
parFdaFitBnd_temp(pp_temp,indFitBnd_temp) = eval_fd(vecfdBnd_temp,pp_temp) ;

% write data
PAR_nadLogRegDkSaFitBnd = parFdaFitBnd_temp ;
par_fdBndCoefs(:,indFitBnd_temp) = getcoef(vecfdBnd_temp) ;

% derive PAR data to get Kd
% compute fd derivative of lum_fd object
KdBndFD_temp = deriv_fd(vecfdBnd_temp) ;
Kd_Bnd(pp_temp,indFitBnd_temp) = eval_fd(KdBndFD_temp,pp_temp) ;

% recover coefficients of computed functional object
% reconstitute array with total dataset size
KdFDcoefs(:,indFitBnd_temp) = getcoef(KdBndFD_temp) ;



%% compute non-log arrays
PAR_nadNonLogRegDkSaFitAll = exp(PAR_nadLogRegDkSaFitAll) ;
PAR_nadNonLogRegDkSaFitBnd = exp(PAR_nadLogRegDkSaFitBnd) ;


%% clear temp data
clear ans *_temp

%% END
