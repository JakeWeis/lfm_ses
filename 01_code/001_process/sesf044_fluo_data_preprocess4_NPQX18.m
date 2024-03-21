%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% compute NPQ non-photochemical quenching correction
% with Xing et al. 2018 method
% 2.1.3 NPQ correction methods based on MLD and light-threshold depth
% WARNING: not expected to perform well in shallow-mixing waters
% NPQ flags:
% 0 OK
% 1 NPQ LAYER MIGHT BE THINNER
% 2 SHALLOW MIXING CONDITIONS
% 7 NIGHT PROFILE
% 8 BOTTOM OF NPQ LAYER NOT REACHED
% 9 NO LIGHT DATA
% method:
% calculate min (PAR15, MLDPhy)
% propagate corresponding max of FLUO within "NPQ" layer up to the surface
% outputs:
% FLUO_nadDkNpq = dark+npq corrected FLUO data

%%%%%%%%%% REFERENCES %%%%%%%%%%
% see Xing et al. 2018 for quenching depth PAR15

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
% last modified: 22.04.18
% 
% 20.03.31 update: adding switch case for plot/no plot option (default =
% YES)
% 21.01.12 update: adjusting variables names (CHLA_nad and chlaData)
% 21.04.26 update: taking MLDPhy into account (instead of MLDbio in
% previous version)
% 22.04.18 update: taking smoothed fluo profile instead of raw profile
% (avoid spikes to be extrapolated to surface leading to high
% overestimates)
% 
% -----------------------------------------------------------------------

disp(strcat('FLUO data pre-processing step5: NPQ / platform:',...
    platform_metadata.platform_name,'.....'))

%% plot MLD and PAR15

switch doPlots
    case 'NO'
        
    case 'YES'

        figure (4201)
        plot(genData.MLDphy(idx041_chlaNonNan),'-k')
        xlabel(strcat('FLUO profile nº (total=',int2str(np041_chlaNonNan),')'))
        ylabel('MLD (m) / PAR15 (m)')
        title(strcat('MLDPhy vs PAR15',' - platform:',platform_metadata.platform_name))
        hold on
        plot(parData.quenchDepth(idx041_chlaNonNan),'-b')
        % highlight where MLD is shallower than PAR15 (~shallow mixing
        % conditions)
        idx_temp = parData.quenchDepth(idx041_chlaNonNan) >...
            genData.MLDphy(idx041_chlaNonNan) ;
        PAR15_values_temp = parData.quenchDepth(idx041_chlaNonNan) ;
        PAR15_values_temp(~idx_temp) = NaN ;
        plot(PAR15_values_temp,'-or')
        legend('MLD','PAR15',...
            strcat('PAR15>MLD (~shallow mixing conditions):',32,int2str(100 * nnz(idx_temp) /...
            nnz(idx041_chlaNonNan)),'%'),...
            'Location','southeast')
        set(gca,'YDir','reverse')
        hold off
        
end



%% apply NPQ correction to FLUO profiles

% choose fluo array to be considered for fluo correction (value to be
% picked up for up-to-the-surface extrapolation)

%%%% moving median of FLUO_nadRegDk array
% width of moving median window
movmedianWindow_temp = 5 ;
% mask of NaN values
maskNan_temp = isnan(FLUO_nadReg) ;
% moving median con FLUO_nadRegDk array
fluoToComputeForNpq_temp = movmedian(FLUO_nadRegDk,movmedianWindow_temp,...
    'omitnan') ;
fluoToComputeForNpq_temp(maskNan_temp) = NaN ;

% initialize FLUO_nadRegDkNpq array
FLUO_nadRegDkNpq = FLUO_nadRegDk ;
% compute min of MLD and PAR15
CORRDEPTH_temp = min(genData.MLDphy,parData.quenchDepth) ;

for ii_temp = find(idx041_chlaNonNan & idx031_lightNonNan).'
    % get index indC_temp corresponding to min(MLD,PAR15)
    indC_temp = find(PRES(:,ii_temp) >= CORRDEPTH_temp(ii_temp),1,'first') ;
    % get index indX12 of zX12
    % corresponding to max of FCHLA within the "NPQ" layer (see Xing et al. 2018)
    [FCHLAzX12_temp,indX12_temp] = max(fluoToComputeForNpq_temp(1:indC_temp,ii_temp)) ;
    % correct FLUO_nadDk
    FLUO_nadRegDkNpq(1:indX12_temp,ii_temp) = FCHLAzX12_temp ;
    % rewrtite effective correcting depth
    if isempty(indX12_temp)
    else
        CORRDEPTH_temp(ii_temp) = indX12_temp ;
    end
end


%% setup a series of flags according to NPQ correction quality control

% initialize flags
flags_temp = nan(np_tot,1) ;

% NPQ LAYER MIGHT BE THINNER
flags_temp(exp(parData.subsurVal) < 15) = 1 ;
% SHALLOW MIXING CONDITIONS
flags_temp(genData.MLDphy(ii_temp) <= parData.quenchDepth) = 2 ;
% NIGHT PROFILE
flags_temp(idx035_lightDay == 0) = 7 ;
% BOTTOM OF NPQ LAYER NOT REACHED
flags_temp(isnan(parData.quenchDepth)) = 8 ;
% NO LIGHT DATA
flags_temp(idx031_lightNonNan == 0) = 9 ;
% OK
flags_temp(isnan(flags_temp)) = 0 ;


%% write data in fluoData table
fluoData.NpqCorrDepth = CORRDEPTH_temp ;
fluoData.NpqFlag = flags_temp ;

%% plot raw (FLUO_nadJUSTED) vs dark+NPQ smoothed corrected profiles


%%%%%%%%%%%%%%%%%% optional plot %%%%%%%%%%%%%%%%%%
switch doPlots
    case 'NO'
        
    case 'YES'

        for ii_temp = find(idx041_chlaNonNan).'
%             pause(1)

            figure (4202)
            % raw FLUO (FLUO_nadJUSTED)
            plot(FLUO_nadReg(:,ii_temp),PRES(:,ii_temp),'-g',...
                'DisplayName','raw FLUO')
            hold on
            % dark+NPQ FLUO
            plot(FLUO_nadRegDkNpq(:,ii_temp),PRES(:,ii_temp),'-k','LineWidth',1,...
                'DisplayName','corrected FLUO (dark noise, NPQ)')
            % edit axes/labels/title/legend
            ax_temp = gca ;
            ax_temp.YDir = 'reverse' ;
            ylim([0 220])
            xlabel('FLUO (mg.m-3)')
            ylabel('depth (m)')
            title({strcat('tag:',platform_metadata.platform_name,...
                    ' - profile nº',int2str(ii_temp),' /',int2str(np_tot)),...
                    strcat('FLUO profile nº',int2str(nnz(idx041_chlaNonNan(1:ii_temp))),...
                    ' /',int2str(np041_chlaNonNan))})
            legend('Location','SouthEast')
            hold off
        end       
end


%% clear temp data
clear ans *_temp

%% END