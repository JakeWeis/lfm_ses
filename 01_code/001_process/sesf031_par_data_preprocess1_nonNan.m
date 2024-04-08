%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% preprocess PAR (photosynthetically available radiation) data
% create parData table to sumarize PAR data

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
% last modified: 22.04.18

% 20.03.31 update: convert lightData_ad matrix into table with col names
% 20.03.31 update: adding switch case for plot/no plot option (default =
% YES)
% 20.04.27 update: compute solar altitude and hour of the day even in light
% NaN profiles
% 20.08.26 update: adding Zeu euphotic depth (Morel and Berthon, 1989)
% 20.09.23 updates:
% removing computation of solar altitude and hour of day (computed in
% genData)
% modifying computation of first/last non NaN values of profile (using
% find_ndim)
% modifying computation of subsur light value subsurVal
% 21.01.12 updates:
% rename lightData_ad variable in lightData
% write metadata (number of non NaN profiles
% 21.02.01 update: modifying computaion of subsur light value subsurVal
% using sub2ind
% 21.04.27 update: adding saturationDepth variable in lightData array
% 21.08.09 update: add switch case for sealtag/float (log/exp data as input)
% 22.04.06 updates:
% renaming lightData into parData
% removing solAlt and HofDay from parData table
% 22.04.18 update: adding meanKd (mean attenuation coefficient) in parData table
% 
% -----------------------------------------------------------------------

disp(strcat('PAR data pre-processing step1: nonNan / platform:',...
    platform_metadata.platform_name,'.....'))

%% create parData table with profile info

% short names
names_temp = cell(1,13) ;
names_temp{1} = 'profileNB' ; % profile number
names_temp{2} = 'firstNonNan' ; % depth of first non NaN value of light
names_temp{3} = 'lastNonNan' ; % depth of last non NaN value of light
names_temp{4} = 'subsurVal' ; % subsurface light value (first non NaN value of light)
names_temp{5} = 'darkDepth' ; % starting depth of dark value (organelli 2016 method)
names_temp{6} = 'darkVal' ; % dark value (organelli 2016 method)
names_temp{7} = 'quenchDepth' ; % depth of quenching threshold 15 umol.m-2.s-1 (xing et al. 2018)
names_temp{8} = 'Zeu' ; % euphotic depth (Morel and Berthon, 1989)
names_temp{9} = 'saturationDepth' ; % signal is valid starting from saturationDepth
names_temp{10} = 'meanKd' ; % mean attenuation coefficient Kd (calculated where LIGHT is non DARK non SAT)
names_temp{11} = 'attSlopTOT' ; % slope of attenuation (calculated where LIGHT is non DARK)
names_temp{12} = 'attSlopPART' ; % slope of attenuation calculated on [min_topIntegBound:min_botIntegBound] interval
names_temp{13} = 'predIntChla' ; % prediction of integrated CHLA from light attenuation slope linear regression

parData = horzcat([1:np_tot].', zeros (np_tot,numel(names_temp)-1)) ;
parData = array2table(parData,'VariableNames',{names_temp{1:end}}) ;




%% compute parData data
switch platform_metadata.platform_type   
    case 'sealtag'
        larray_temp = PAR_nadLogReg ;
    case 'float'
        larray_temp = PAR_nadNonLogReg ;
end

% first/last nonNaN value of each profile (based on PAR_nadLogReg array)
parData.firstNonNan = find_ndim(~isnan(larray_temp),1,'first').' ;
parData.firstNonNan(parData.firstNonNan == 0) = NaN ;
parData.lastNonNan = find_ndim(~isnan(larray_temp),1,'last').' ;
parData.lastNonNan(parData.lastNonNan == 0) = NaN ;

% subsurface value
% subsurIndexes_temp = sub2ind(size(larray_temp),...
%     parData.firstNonNan,...
%     transpose(1:platform_metadata.np)) ;
% subsurValues_temp = nan(platform_metadata.np,1) ;
% subsurValues_temp(~isnan(subsurIndexes_temp)) =...
%     larray_temp(subsurIndexes_temp(~isnan(subsurIndexes_temp))) ;
subsurValues_temp = NaN(platform_metadata.np,1);
finiteIndices = find(isfinite(parData.firstNonNan));
subsurValues_temp(finiteIndices,1) = arrayfun(@(x) larray_temp(parData.firstNonNan(x), x), finiteIndices);

switch platform_metadata.platform_type   
    case 'sealtag'
        parData.subsurVal = subsurValues_temp ;
    case 'float'
        subsurValues_temp(subsurValues_temp <= 0) = NaN ;
        subsurValues_temp = log(subsurValues_temp) ;
end
parData.subsurVal = subsurValues_temp ;


%% create logical indexing for nonNaN profiles

% non NaN LIGHT profiles 
idx031_lightNonNan = ~isnan(parData.firstNonNan(:)) ;
np031_lightNonNan = nnz(idx031_lightNonNan) ;
disp(strcat('dataset:',...
    int2str(np031_lightNonNan),...
    '/',...
    int2str(platform_metadata.np),...
    ' non NaN LIGHT profiles'))



%% proceed to visual check to verify consistency of selected profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: CAN BE FASTIDUOUS (LARGE NUMBER OF PROFILES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch doPlots
    case 'NO'
        
    case 'YES'
        for ii_temp = 1:np_tot
            if idx031_lightNonNan(ii_temp) == 0
            else
                figure (311)
                plot(PAR_nadLogReg(:,ii_temp),PRES(:,ii_temp),'-k')
                % can choose to plot all profiles on same figure
                % hold on
                xlabel('LIGHT (ln(umol.m-2.s-1))')
                ylabel('PRES (decibar)')
                title({strcat('tag:',platform_metadata.platform_name,...
                    ' - profile nº',int2str(ii_temp),' /',int2str(np_tot)),...
                    strcat('LIGHT profile nº',int2str(nnz(idx031_lightNonNan(1:ii_temp))),...
                    ' /',int2str(np031_lightNonNan))})
                %pause(0.5)
            end
        end

end



%% write metadata
platform_metadata.np031_lightNonNan = nnz(idx031_lightNonNan) ;

%% clear temp data
clear ans *_temp


%% END
