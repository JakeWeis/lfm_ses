%% script metadata

% ----------------------------------------------------------------------
%%%%%%%%%% SCRIPT ACTIONS %%%%%%%%%%
% preprocess PAR data
% detect saturated values of PAR (where dPAR/dPRES(surface)  >= 0)
% (surface values, high light conditions)
% erase profile values where saturated


%%%%%%%%%% DEPENDANT SCRIPTS (previously created variables) %%%%%%%%%%
% sesf000_define_platform
% sesf00x
% sesf01x
% sesf03x

%%%%%%%%%% MATLAB VERSION %%%%%%%%%%
disp('***MATLAB version information***')
disp(strcat('script: MATLAB version 9.10.0.1669831 (R2021a) Update 2'))
disp(['current MATLAB version:' version])

%%%%%%%%%% AUTHOR / LAST MODIFIED %%%%%%%%%%
% L. Le Ster (lls)
% last modified: 22.04.18
% 
% 22.04.06 update: renaming LIGHT into PAR (L_adxxx becomes PAR_nadxxx)
% 22.04.15 update: condition for saturation is not only PAR profile being
% constant or increasing with depth (derivative >= 0) but decreasing less
% than mean decreasing rate (derivative < mean derivative) - derivative is
% computed on log profile i.e. we take mean Kd into account
%
% -----------------------------------------------------------------------

disp(strcat('PAR data pre-processing step3: sensor saturation / platform:',...
    platform_metadata.platform_name,'.....'))

%% detect saturation depth + correct PAR profile by removing saturated portions

% initialize data
PAR_nadNonLogRegDkSa = PAR_nadNonLogRegDk ;
satDepth_temp = ones(np_tot,1) ; % (default value = 1 i.e. profile is
% valid from surface)

for iProfile_temp = find(idx031_lightNonNan).'
    % compute derivative of PAR_nadNonLogRegDk profile
    diffPARlog_temp = diff(PAR_nadLogRegDk(:,iProfile_temp)) ;
    % mean of local derivative of PAR_nadNonLogRegDk profile
    meandiffPARlog_temp = mean(diffPARlog_temp,'omitnan') ;
    % check saturation depth only if profile starts with
    % a portion with descreasing intensity lower than mean derivative (i.e.
    % constant/increasing/low decreasing portion)
    if diffPARlog_temp(parData.firstNonNan(iProfile_temp)) < meandiffPARlog_temp
    else
        % initialize satDepth_temp
        satTestDepth_temp = parData.firstNonNan(iProfile_temp) + 1 ;
        while diffPARlog_temp(satTestDepth_temp) >= meandiffPARlog_temp
            satTestDepth_temp = satTestDepth_temp + 1 ;
        end
        % signal is valid from satDepth_temp
        satDepth_temp(iProfile_temp) = satTestDepth_temp ;
        % rewrite PAR_nadNonLogRegDkSa to remove saturated values
        PAR_nadNonLogRegDkSa(1:satTestDepth_temp - 1,iProfile_temp) = NaN ;
    end
end

%% write data in parData table

parData.saturationDepth = satDepth_temp ;


%% compute natural log of dark+saturation corrected PAR data
PAR_nadLogRegDkSa = log(PAR_nadNonLogRegDkSa) ;


%% clear temp data

clear ans *_temp

%% END
    