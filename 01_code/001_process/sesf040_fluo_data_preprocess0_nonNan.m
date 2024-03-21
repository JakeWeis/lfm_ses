%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% create fluoData array with FLUO summarized data for each profile
% remove test / inconsitent profiles (indexes set by user in rem_temp vector)
% remove bottom part of profiles (i.e. depth > maxCHLA_depth)

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
% last modified: 22.03.29

% 20.02.11 update:
% adding 1 line in chlaData_ad array for MLD (physical)
% new chlaData_ad array size = 9 x np_tot
% 20.02.28 update:
% adding 1 line in chlaData_ad array for integrated CHLA (satellite x MLD physical)
% new chlaData_ad array size = 10 x np_tot
% 20.03.30 update : convert chlaData_ad array into table with col names
% 20.03.31 update: adding switch case for plot/no plot option (default =
% YES)
% 20.09.23 update: modifying computation of chlaData_ad variables
% 21.01.12 updates:
% changing chlaData_ad variable name to chlaData
% adding subsurVal and surfChla distinct variables in chlaData table
% changing idx021 to idx041 and np021 to np041 (naming consistency)
% rewrite CHLA_nad array to rewrite with NaN all inconsistent profiles
% 21.02.01 update: modifying computation of subsur chla value using sub2ind
% 21.03.31 update: modifying ylim in plot to focus on [0 maxCHLA_depth]
% 21.04.26 updates:
% adapting chlaData table variables names to new processing steps
% moving first non Nan / last non Nan computations to sesf041 script
% 22.02.03 update: reversing yaxis of chla plot
% 22.03.29 update: chlaData table becomes fluoData table
%
% -----------------------------------------------------------------------

disp(strcat('FLUO data pre-processing step0: nonNan / platform:',...
    platform_metadata.platform_name,'.....'))


%% create fluoData table

% variables names
names_temp = {} ;
names_temp(1) = {'profileNB'} ; % profile number
% profile metrics (raw data)
names_temp(2) = {'firstNonNan'} ; % depth of first non NaN value of fluo
names_temp(3) = {'lastNonNan'} ; % depth of last non NaN value of fluo
% profile correction
names_temp(4) = {'darkValProfile'} ; % dark value (median of 190 m-200 m values)
names_temp(5) = {'darkValCorr'} ; % dark value for correction (movmedian + linear fit)
names_temp(6) = {'darkDepth'} ; % depth of first value belonging to dark signal
names_temp(7) = {'RsdMixing'} ; % RSD relative standard deviation within mixing layer for quality control
names_temp(8) = {'NpqCorrDepth'} ; % depth from which NPQ correction is applied (see sesf045 script)
names_temp(9) = {'NpqFlag'} ; % NPQ correction flag (see sesf045 script)
% profile metrics (post dark+npq correction)
names_temp(10) = {'maxVal'} ; % fluo max value
names_temp(11) = {'maxValDepth'} ; % depth of fluo max value
names_temp(12) = {'MLDbio'} ; % MLDbio (Lacour et al. 2017)
names_temp(13) = {'surf'} ; % surface value of fluo
names_temp(14) = {'integ'} ; % integrated fluo (trapezoidal integration)
names_temp(15) = {'smoothRelResid'} ; % smoothing residuals (%) 

fluoData = horzcat([1:np_tot].', nan(np_tot,numel(names_temp)-1)) ;
fluoData = array2table(fluoData,'VariableNames',{names_temp{1:end}}) ;


%% rewrite FLUO data to remove values below max sampled depth
% (maxCHLA_depth is set by platform operator before deployment)

FLUO_nadReg(PRES(:) > maxCHLA_depth) = NaN ;


%% detect empty profiles

% 1º) nan profiles
idxEmptyNan_temp = find_ndim(~isnan(FLUO_nadReg),1) == 0 ;
% 2º) according to abs(max-min) value of profile (range)
% ses empty FLUO profiles are filled with zeros or constant values
idxEmptyZeros_temp = range(FLUO_nadReg,1) == 0 ;

%% detect profiles defined on narrow interval
% (too narrow for functional fit)
numElPerProfile_temp = arrayfun(@(a) nnz(~isnan(FLUO_nadReg(:,a))),1:platform_metadata.np) ;
idxTooNarrowForfFit_temp = numElPerProfile_temp <= nbaz ;

%% create logical indexing for nonNaN profiles

% non NaN FLUO profiles 
idx041_chlaNonNan = (~idxEmptyNan_temp & ~idxEmptyZeros_temp &...
    ~idxTooNarrowForfFit_temp ).' ;
np041_chlaNonNan = nnz(idx041_chlaNonNan) ;
disp(strcat('dataset:',...
    int2str(np041_chlaNonNan),...
    '/',...
    int2str(platform_metadata.np),...
    ' non NaN raw FLUO profiles'))


        
%% proceed to visual check to verify consistency of selected profiles
switch doPlots
    case 'NO'
        
    case 'YES'
        lineCol_temp = cell(1,np_tot) ;
        lineCol_temp(idx035_lightDay) = {'-g'} ;
        lineCol_temp(~idx035_lightDay) = {'-b'} ;
        for ii_temp = find(idx041_chlaNonNan).'
            figure (2101)
            plot(FLUO_nadReg(:,ii_temp),PRES(:,ii_temp),...
                char(lineCol_temp(ii_temp)))
            xlabel('FLUO (mg.m-3)')
            ylabel('PRES (decibar)')
            title({strcat('tag:',platform_metadata.platform_name,...
                ' - profile nº',int2str(ii_temp),' /',int2str(np_tot)),...
                strcat('FLUO profile nº',int2str(nnz(idx041_chlaNonNan(1:ii_temp))),...
                ' /',int2str(np041_chlaNonNan))})
            ylim([0 maxCHLA_depth])
            xlim([0 5])
            set(gca,'YDir','reverse')
            grid on
            grid minor
            pause(0.5)
        end

end

%% display information about number of non NaN profiles

np041_chlaNonNan = nnz(idx041_chlaNonNan);
disp(strcat('dataset:',...
    int2str(np041_chlaNonNan),...
    '/',...
    int2str(platform_metadata.np),...
    ' non NaN FLUO profiles (test profiles removed)'))

%% rewrite data
FLUO_nadReg(:,~idx041_chlaNonNan) = NaN ;


%% write metadata
platform_metadata.np041_chlaNonNan = np041_chlaNonNan ;


%% clear temp data
clear *_temp


%% END
