%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% compute RSD relative standard deviation of fluo signal within
% mixing layer and below the quenching depth (quenching depth<.<MLD)
% for quality control purposes
% also known as coefficient of variation (CV)
% e.g. all winter profiles with RSD > 20% removed in Lacour et al. 2017

%%%%%%%%%% REFERENCES %%%%%%%%%%
% RSD: Lacour et al. 2017 additional data
% see Xing et al. 2018 for quenching depth

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
% J. Weis (jw)
% last modified: 22.04.06
% 
% 21.01.12 update: adjusting variables names (CHLA_nad and chlaData)
% 21.04.26 update: consider MLDPhy (instead of MLDbio in previous version)
% 
% -----------------------------------------------------------------------

disp(strcat('FLUO data pre-processing step8: RSD / platform:',...
    platform_metadata.platform_name,'.....'))

%% compute RSD relative standard deviation

rsd_temp = nan(np_tot,1) ;

idx_temp = idx041_chlaNonNan & idx031_lightNonNan & ~isnan(genData.MLDphy) ;
for ii_temp = find(idx_temp).'
    if genData.MLDphy(ii_temp) < parData.quenchDepth(ii_temp) ||...
            isnan(parData.quenchDepth(ii_temp))
        topInd_temp = 1 ;
    else
        topInd_temp = parData.quenchDepth(ii_temp) ;
    end
    botInd_temp = genData.MLDphy(ii_temp) ;
    % standard deviation of CHLA in topInd_temp:botInd_temp interval
    stdFluo_temp = std(FLUO_nadRegDk(topInd_temp:botInd_temp,ii_temp)) ;
    % mean of CHLA in topInd_temp:botInd_temp interval
    meanFluo_temp = mean(FLUO_nadRegDk(topInd_temp:botInd_temp,ii_temp)) ;
    if meanFluo_temp == 0
    else
        % computation of RSD relative standard deviation
        rsd_temp(ii_temp) = stdFluo_temp / meanFluo_temp ;
    end
end

% write nan elsewhere
fluoData.RsdMixing = rsd_temp ;

    
%% clear temp data

clear ans *_temp

%% END