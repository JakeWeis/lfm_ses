%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% compute Mixed Layed Depth (physical, using TEMP and SAL)
% based on threshold criteria deltad = 0.03 kg.m-3

%%%%%%%%%% REFERENCES %%%%%%%%%%
% de Boyer Montégut et al. 2004 doi:10.1029/2004JC002378.

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
% last modified: 22.04.18
% 
% 22.04.06 update: compute MLD003 on DENS_fit instead of DENS raw (DENS raw
% is computed on non-ADJUSTED fields of TEMP and SAL and is often noisy)
% 22.04.18 update: no more DENS_fit computation => MLD003 is computed on
% DENS array (DENS is computed with TEMP_ADJUSTED and PSAL_ADJUSTED fields)
% 
% -----------------------------------------------------------------------

disp(strcat('computing MLD for platform:',...
    platform_metadata.platform_name,'.....'))

%% compute MLD using method: threshold method (0.03 kg.m-3)

% DENS array to be used for MLD003 computation
densArrayforMLDComputation_temp = DENS ;

% reference values
refDepth_temp = 10 ;
minProfileStartDepth_temp = 20 ;
densThreshold_temp = 0.03 ;

% filter valuable data
idxMinProfileStart_temp = genData.firstNonNan <= minProfileStartDepth_temp ;

% index of first value below refDepth depth (default = 10 m)
refDepthIdxVec_temp = find_ndim(PRES >= refDepth_temp & ~isnan(TEMP),1,'first').' ;
refDepthIdxVec_temp(refDepthIdxVec_temp == 0) = NaN ;
% reference depth value
refDepthIndexes_temp = sub2ind(size(PRES),...
    refDepthIdxVec_temp,...
    transpose(1:np_tot)) ;
refDepthValues_temp = nan(np_tot,1) ;
refDepthValues_temp(~isnan(refDepthIndexes_temp)) =...
    PRES(refDepthIndexes_temp(~isnan(refDepthIndexes_temp))) ;
refDepthValues_temp(~idxMinProfileStart_temp) = NaN ;

% reference density value
refDensValues_temp = nan(np_tot,1) ;
refDensValues_temp(~isnan(refDepthIndexes_temp)) =...
    densArrayforMLDComputation_temp(refDepthIndexes_temp(~isnan(refDepthIndexes_temp))) ;

% compute sigma0
% (potential density referenced to the surface)
sigma0_temp = densArrayforMLDComputation_temp - repmat(refDensValues_temp.',max_depth,1) ;

% compute MLD
% mld index
mldIndex_temp = find_ndim(sigma0_temp > densThreshold_temp &...
    PRES > refDepth_temp, 1, 'first').' ;
mldIndex_temp(mldIndex_temp == 0) = NaN ;
% mld depth value
mldIndexes_temp = sub2ind(size(PRES),...
    mldIndex_temp,...
    transpose(1:np_tot)) ;
mldValues_temp = nan(np_tot,1) ;
mldValues_temp(~isnan(mldIndexes_temp)) =...
    PRES(mldIndexes_temp(~isnan(mldIndexes_temp))) ;


genData.MLD003 = mldValues_temp ;

%% plot MLD series for entire platform transect

switch doPlots
    case 'NO'
        
    case 'YES'
        figure (181)
        plot1_temp = plot(genData.MLD003,'.','MarkerSize',10,'Color', [0 0 1]*150/255,...
            'DisplayName','MLD003') ;
        xlabel('profile nº')
        ylabel('MLD (m)')
        hold on
        plot2_temp = plot(movmedian(genData.MLD003,180,'omitnan'),'-c','LineWidth',2,...
            'DisplayName','moving median') ;
        title(strcat('series of MLD',' - platform:',platform_metadata.platform_name,...
            ' (',int2str(np_tot),' DENSITY profiles)'))
        legend()
        grid on
        grid minor
        hold off

end



%% clear temp variables
clear *_temp


%% END
