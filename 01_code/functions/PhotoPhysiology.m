function [Data,ProfileInfo] = PhotoPhysiology(Data,ProfileInfo,defaultPars)
%% CMD message: start
if isempty(getCurrentTask)
    fprintf('Calculating <strong>photo-physiological metrics</strong>...');
end

%% Create parData info table
var_names = {...
    'Profile', ...              % profile number
    'noData_Fluo',...           % no data in profile
    'noData_Light',...          % no data in profile
    'depth_upper',...           % upper limit of depth range based on which metrics are calculated
    'depth_lower',...           % lower limit of depth range based on which metrics are calculated
    'fluo_mean',...             % mean fluorescence
    'fluo_SD',...               % SD fluorescence
    'fluo_integr',...           % integrated fluorescence
    'Kd_bspline',...            % Kd --> mean of depth resolved Kd values calculated in processRadiometry_fit (bspline fitted derivative PAR function)
    'Kd_bspline_SD',...         % Kd --> SD of depth resolved Kd values calculated in processRadiometry_fit (bspline fitted derivative PAR function)
    'Kd_linfit',...             % Kd --> ln(PAR) vs. depth regression slope
    'Kd_linfit_SD',...          % Kd --> ln(PAR) vs. depth standard error
    'Kd_linfit_pval',...        % Kd --> ln(PAR) vs. depth p-value
    'Kd_linfit_R2',...          % Kd --> ln(PAR) vs. depth R2
    'Kd490_bspline',...         % Kd490 from DOWN_IRRADIANCE490
    'Kd490_bspline_SD',...
    'Kd490_linfit',...
    'Kd490_linfit_SD',...
    'Kd490_linfit_pval',...
    'Kd490_linfit_R2',...
    'CHL_Kd490_bspline',...     % Chl from Kd490
    'CHL_Kd490_bspline_SD',...
    'CHL_Kd490_linfit',...
    'CHL_Kd490_linfit_SD',...
    'SlopeFactor_bspline',...   % Slope factor
    'SlopeFactor_bspline_SD',...
    'SlopeFactor_linfit',...
    'SlopeFactor_linfit_SD',...
    'NPQ_alpha',...             % NPQ metrics
    'NPQ_slope',...
    'NPQ_max'}';
ProfileInfo.PhotoPhys = array2table(NaN(Data.MetaData.nProfs, numel(var_names)),'VariableNames',var_names);

%% Populate initial information
% Profile #
ProfileInfo.PhotoPhys.Profile = (1 : Data.MetaData.nProfs)';

% Find profiles without usable light or fluorescence data (all NaNs)
% Light variable name
if isfield(Data.Processed,'LIGHT')
    light_varname = 'LIGHT';
elseif isfield(Data.Processed,'DOWNWELLING_PAR')
    light_varname = 'DOWNWELLING_PAR';
elseif isfield(Data.Processed,'DOWN_IRRADIANCE490')
    light_varname = 'DOWN_IRRADIANCE490';
end
ProfileInfo.PhotoPhys.noData_Fluo = all(isnan(Data.Processed.FLUO.Reg))' | range(Data.Processed.FLUO.Reg)' == 0;
ProfileInfo.PhotoPhys.noData_Light = all(isnan(Data.Processed.(light_varname).log.Reg))' | range(Data.Processed.(light_varname).log.Reg)' == 0;

% depth range based on which metrics are calculated
z_range = defaultPars.PhotoPhys.depthrange;
ProfileInfo.PhotoPhys.depth_upper(:) = z_range(1);
ProfileInfo.PhotoPhys.depth_lower(:) = z_range(end);

%% Fluorescence (mean, integrated)
fluo_range = max(Data.Processed.FLUO.RegDrkNPQFitAll(z_range,:),0);
fluo_range(:,all(fluo_range==0)) = NaN;
ProfileInfo.PhotoPhys.fluo_mean = mean(fluo_range,1,'omitnan')';
ProfileInfo.PhotoPhys.fluo_SD = std(fluo_range,0,1,'omitnan')';
for iProf = 1 : Data.MetaData.nProfs
    fluo_iProf = fluo_range(:,iProf);
    ProfileInfo.PhotoPhys.fluo_integr(iProf,1) = trapz(z_range(isfinite(fluo_iProf)),fluo_iProf(isfinite(fluo_iProf)),1);
end
ProfileInfo.PhotoPhys.fluo_mean(ProfileInfo.PhotoPhys.fluo_mean <= 0) = NaN;
ProfileInfo.PhotoPhys.fluo_integr(ProfileInfo.PhotoPhys.fluo_integr <= 0) = NaN;

%% Attenuation from LIGHT (seal tags) or DOWNWELLING_PAR (float)
if isfield(Data.Processed,'LIGHT') || isfield(Data.Processed,'DOWNWELLING_PAR')
    % Kd --> mean of depth resolved Kd values calculated in processRadiometry_fit (bspline fitted derivative PAR function)
    Kd_bspline_range = Data.Processed.(light_varname).Kd.FitAll(z_range,:);
    ProfileInfo.PhotoPhys.Kd_bspline = mean(Kd_bspline_range,'omitnan')';
    ProfileInfo.PhotoPhys.Kd_bspline_SD = std(Kd_bspline_range,0,'omitnan')';

    % Kd --> linear regression slope of ln(PAR) vs. depth
    % lnlight_range = Data.Processed.(light_varname).log.RegDrkSatFitAll(z_range,:);
    lnlight_range = Data.Processed.(light_varname).log.RegDrkSat(z_range,:);
    for iProf = 1 : Data.MetaData.nProfs
        if numel(find(isfinite(lnlight_range(:,iProf))))>2
            lfit = fitlm(z_range,lnlight_range(:,iProf));
            ProfileInfo.PhotoPhys.Kd_linfit(iProf,1) = -lfit.Coefficients.Estimate(2);
            ProfileInfo.PhotoPhys.Kd_linfit_SD(iProf,1) = -lfit.Coefficients.SE(2);
            ProfileInfo.PhotoPhys.Kd_linfit_pval(iProf,1) = lfit.Coefficients.pValue(2);
            ProfileInfo.PhotoPhys.Kd_linfit_R2(iProf,1) = lfit.Rsquared.Ordinary;
        end
    end
end

%% Attenuation from DOWN_IRRADIANCE490
if isfield(Data.Processed,'DOWN_IRRADIANCE490')
    % Kd --> mean of depth resolved Kd values calculated in processRadiometry_fit (bspline fitted derivative PAR function)
    Kd490_bspline_range = Data.Processed.DOWN_IRRADIANCE.Kd.FitAll(z_range,:);
    ProfileInfo.PhotoPhys.Kd490_bspline = mean(Kd490_bspline_range,'omitnan')';
    ProfileInfo.PhotoPhys.Kd490_bspline_SD = std(Kd490_bspline_range,0,'omitnan')';

    % Kd --> linear regression slope of ln(PAR) vs. depth
    % lnlight_range = Data.Processed.DOWN_IRRADIANCE.log.RegDrkSatFitAll(z_range,:);
    lnlight_range = Data.Processed.DOWN_IRRADIANCE.log.RegDrkSat(z_range,:);
    for iProf = 1 : Data.MetaData.nProfs
        if numel(find(isfinite(lnlight_range(:,iProf))))>2
            lfit = fitlm(z_range,lnlight_range(:,iProf));
            ProfileInfo.PhotoPhys.Kd490_linfit(iProf,1) = -lfit.Coefficients.Estimate(2);
            ProfileInfo.PhotoPhys.Kd490_linfit_SD(iProf,1) = -lfit.Coefficients.SE(2);
            ProfileInfo.PhotoPhys.Kd490_linfit_pval(iProf,1) = lfit.Coefficients.pValue(2);
            ProfileInfo.PhotoPhys.Kd490_linfit_R2(iProf,1) = lfit.Rsquared.Ordinary;
        end
    end
else
    % Calulate Kd_490 from Kd_PAR (float) or Kd_Light (seal tag)
    [ProfileInfo.PhotoPhys.Kd490_bspline,ProfileInfo.PhotoPhys.Kd490_bspline_SD] = Kd_to_Kd490(ProfileInfo.PhotoPhys.Kd_bspline,ProfileInfo.PhotoPhys.Kd_bspline_SD,1);
    [ProfileInfo.PhotoPhys.Kd490_linfit,ProfileInfo.PhotoPhys.Kd490_linfit_SD] = Kd_to_Kd490(ProfileInfo.PhotoPhys.Kd_linfit,ProfileInfo.PhotoPhys.Kd_linfit_SD,1);
end

% Set KD490 minimum to 0.0166 (attenuation in pure water)
ProfileInfo.PhotoPhys.Kd490_bspline(ProfileInfo.PhotoPhys.Kd490_bspline<.0166) = .0166;
ProfileInfo.PhotoPhys.Kd490_linfit(ProfileInfo.PhotoPhys.Kd490_linfit<.0166) = .0166;

%% Chl from Kd490 (±SD)
% following Morel et al. (2007), equ. 8, https://doi.org/10.1016/j.rse.2007.03.012
ProfileInfo.PhotoPhys.CHL_Kd490_bspline = ((ProfileInfo.PhotoPhys.Kd490_bspline - 0.0166)/0.0773).^(1/0.6715);
ProfileInfo.PhotoPhys.CHL_Kd490_bspline_SD = (ProfileInfo.PhotoPhys.CHL_Kd490_bspline*(1/0.6715).*(ProfileInfo.PhotoPhys.Kd490_bspline_SD/0.0773))./((ProfileInfo.PhotoPhys.Kd490_bspline - 0.0166)/0.0773);
ProfileInfo.PhotoPhys.CHL_Kd490_linfit = ((ProfileInfo.PhotoPhys.Kd490_linfit - 0.0166)/0.0773).^(1/0.6715);
ProfileInfo.PhotoPhys.CHL_Kd490_linfit_SD = (ProfileInfo.PhotoPhys.CHL_Kd490_linfit*(1/0.6715).*(ProfileInfo.PhotoPhys.Kd490_linfit_SD/0.0773))./((ProfileInfo.PhotoPhys.Kd490_linfit - 0.0166)/0.0773);
% chl0 = ProfileInfo.PhotoPhys.CHL_Kd490_bspline==0;
% ProfileInfo.PhotoPhys.CHL_Kd490_bspline(chl0) = NaN;
% ProfileInfo.PhotoPhys.CHL_Kd490_bspline_SD(chl0) = NaN;
% chl0 = ProfileInfo.PhotoPhys.CHL_Kd490_linfit==0;
% ProfileInfo.PhotoPhys.CHL_Kd490_linfit(chl0) = NaN;
% ProfileInfo.PhotoPhys.CHL_Kd490_linfit_SD(chl0) = NaN;

% Remove CHL_KD490 values that are lower than the minimum FCHL_MEAN value
% lowVals = ProfileInfo.PhotoPhys.CHL_Kd490_bspline<min(ProfileInfo.PhotoPhys.fluo_mean);
% ProfileInfo.PhotoPhys.CHL_Kd490_bspline(lowVals) = NaN;
% ProfileInfo.PhotoPhys.CHL_Kd490_bspline_SD(lowVals) = NaN;
% lowVals = ProfileInfo.PhotoPhys.CHL_Kd490_linfit<min(ProfileInfo.PhotoPhys.fluo_mean);
% ProfileInfo.PhotoPhys.CHL_Kd490_linfit(lowVals) = NaN;
% ProfileInfo.PhotoPhys.CHL_Kd490_linfit_SD(lowVals) = NaN;

%% Slope factor
% following Schallenberg et al. (2022)
ProfileInfo.PhotoPhys.SlopeFactor_bspline = ProfileInfo.PhotoPhys.fluo_mean./ProfileInfo.PhotoPhys.CHL_Kd490_bspline;
ProfileInfo.PhotoPhys.SlopeFactor_bspline_SD = ProfileInfo.PhotoPhys.SlopeFactor_bspline .* sqrt((ProfileInfo.PhotoPhys.fluo_SD./ProfileInfo.PhotoPhys.fluo_mean).^2+(ProfileInfo.PhotoPhys.CHL_Kd490_bspline_SD./ProfileInfo.PhotoPhys.CHL_Kd490_bspline).^2);
SFinf = ProfileInfo.PhotoPhys.SlopeFactor_bspline == abs(inf);
ProfileInfo.PhotoPhys.SlopeFactor_bspline(SFinf) = NaN;
ProfileInfo.PhotoPhys.SlopeFactor_bspline_SD(SFinf) = NaN;
ProfileInfo.PhotoPhys.SlopeFactor_linfit = ProfileInfo.PhotoPhys.fluo_mean./ProfileInfo.PhotoPhys.CHL_Kd490_linfit;
ProfileInfo.PhotoPhys.SlopeFactor_linfit_SD = ProfileInfo.PhotoPhys.SlopeFactor_linfit .* sqrt((ProfileInfo.PhotoPhys.fluo_SD./ProfileInfo.PhotoPhys.fluo_mean).^2+(ProfileInfo.PhotoPhys.CHL_Kd490_linfit_SD./ProfileInfo.PhotoPhys.CHL_Kd490_linfit).^2);
SFinf = ProfileInfo.PhotoPhys.SlopeFactor_linfit == abs(inf);
ProfileInfo.PhotoPhys.SlopeFactor_linfit(SFinf) = NaN;
ProfileInfo.PhotoPhys.SlopeFactor_linfit_SD(SFinf) = NaN;

%% alpha_NPQ and initial NPQ slope
NPQprofs = find(isfinite(ProfileInfo.FLUO.NPQDepth) & ProfileInfo.FLUO.NPQDepth < -10)';
% ct = 0;
for iProf = NPQprofs
    % ct = ct+1;
    light_range = Data.Processed.(light_varname).lin.RegDrkSatFitSurf(:,iProf);
    NPQ = (Data.Processed.FLUO.RegDrkNPQ(:,iProf) - Data.Processed.FLUO.RegDrk(:,iProf)) ./ Data.Processed.FLUO.RegDrk(:,iProf);
    NPQ(NPQ==0) = NaN;
    finiteVals = isfinite(light_range) & isfinite(NPQ);
    if numel(find(finiteVals))>1
        light_range(~finiteVals) = [];
        NPQ(~finiteVals) = [];

        %% NPQ vs. PAR exponential fit
        % Offsetting PAR by PAR at first NPQ depth to force fit through different origin
        PAR_NPQstart = light_range(end);
        PAR_shift = light_range - PAR_NPQstart;
        % Define fit type
        fitType = fittype('NPQ_max*(1-exp(-(PAR*alpha_NPQ)/NPQ_max))', 'independent', 'PAR', 'coefficients', {'alpha_NPQ', 'NPQ_max'});
        NPQ_max = max(NPQ);
        % Set fit options
        fitOptions = fitoptions(fitType);
        % [alpha_NPQ, NPQ_max]
        fitOptions.StartPoint = [0.01, NPQ_max];
        fitOptions.Lower = [0, 0];
        fitOptions.Upper = [5, 20];
        % Fit model to data
        [fitResult_exp, ~] = fit(PAR_shift, NPQ, fitType, fitOptions);

        % Store alpha_NPQ model parameter
        ProfileInfo.PhotoPhys.NPQ_alpha(iProf,1) = fitResult_exp.alpha_NPQ;
        ProfileInfo.PhotoPhys.NPQ_max(iProf,1) = fitResult_exp.NPQ_max;

        %% Initial NPQ vs. PAR linear regression slope
        % PAR limit for linear regression calculation (NPQ flattens at high PAR)
        % NB: float PAR appears to be considerably higher than seal tag PAR, thus different limits apply
        if strcmp(Data.MetaData.platform_type,'sealtag')
            PARlimit = inf; % set infinite to consider all values...?
        elseif strcmp(Data.MetaData.platform_type,'float')
            PARlimit = inf;
        end

        fitResult_lin = fitlm(light_range(light_range<PARlimit),NPQ(light_range<PARlimit),'Intercept',false);
        R = fitResult_lin.Rsquared.Ordinary;
        m = fitResult_lin.Coefficients.Estimate(1);
        if R>0.25
            ProfileInfo.PhotoPhys.NPQ_slope(iProf,1) = m(R==max(R));
        end
    end

    %% R&D PLAYGROUND: Some more experimental NPQ metrics reliant only on F... wish me luck.
    % NPQ = (Data.Processed.FLUO.RegDrkNPQ(:,iProf) - Data.Processed.FLUO.RegDrk(:,iProf)) ./ Data.Processed.FLUO.RegDrk(:,iProf);
    % NPQ(NPQ==0) = NaN;
    % finiteVals = isfinite(NPQ);
    % if numel(find(finiteVals))>1
    %     i_fval = find(finiteVals);
    %     OUT_NPQsurf(ct,1) = NPQ(i_fval(1));
    %     OUT_NPQarea(ct,1) = trapz(find(finiteVals),NPQ(finiteVals));
    %     % OUT_NPQareaPAR(ct,1) = trapz(PAR,NPQ);
    %     OUT_NPQphi(ct,1) = OUT_NPQarea(ct)./range(find(finiteVals));
    %     OUT_NPQangle(ct,1) = OUT_NPQsurf(ct)./sqrt(OUT_NPQsurf(ct).^2+range(find(finiteVals)).^2);
    % end
end

%% CMD message: done
if isempty(getCurrentTask)
    fprintf('\b\b \x2713\n')
end

%% Function to convert Kd to Kd_490
    function [Kd490,s_Kd490] = Kd_to_Kd490(KdPAR,s_KdPAR,equation)
        % Calculating Kd490 from KdPAR (and uncertainties)
        % after Morel et al. (2007), equs. 9 and 9', https://doi.org/10.1016/j.rse.2007.03.012
        % Equations are given with KdPAR as the unknown and need to be rearranged to be used to calculate Kd490

        % 1) Equation:
        % KdPAR = X + Y*Kd490 − Z*Kd490^−1
        % where X = [0.0864,0.0665], Y = [0.884,0.874], Z = [-0.00137,-0.00121]

        % 2) Rearrange into quadratic equation (a*x^2 + b*x + c = 0, where x = Kd490):
        % Y*Kd490^2 + (X-KdPAR)*Kd490 - Z = 0
        % Quadratic equation constants: a = Y, b = X-KdPAR, c = Z

        % 3) Solve for Kd490 (x = (-b±sqrt(b^2-4ac))/2a):
        % Kd490 = (-(X-KDPAR)±sqrt((X-KDPAR)^2+4YZ))/2Y


        % Equation 1:
        a = 0.884;
        s_a = 0;
        b = 0.0864 - KdPAR;
        s_b = s_KdPAR;
        c = -0.00137;
        s_c = 0;
        % Calculate the solutions
        discriminant = b.^2 - 4*a*c;
        s_discriminant = (b.^2.*2.*s_b)./b;
        Kd490_1a = (-b + sqrt(discriminant)) ./ (2*a);
        s_sqrt_discr = (sqrt(discriminant).*0.5.*s_discriminant)./discriminant;
        s_Kd490_1a = sqrt(s_b.^2 + s_sqrt_discr.^2) ./ (2*a);
        % Kd490_1b = (-b - sqrt(discriminant)) ./ (2*a);

        % Equation 2:
        a = 0.874;
        s_a = 0;
        b = 0.0665 - KdPAR;
        s_b = s_KdPAR;
        c = -0.00121;
        s_c = 0;
        % Calculate the solutions
        discriminant = b.^2 - 4*a*c;
        s_discriminant = (b.^2.*2.*s_b)./b;
        Kd490_2a = (-b + sqrt(discriminant)) ./ (2*a);
        s_sqrt_discr = (sqrt(discriminant).*0.5.*s_discriminant)./discriminant;
        s_Kd490_2a = sqrt(s_b.^2 + s_sqrt_discr.^2) ./ (2*a);
        % Kd490_2b = (-b - sqrt(discriminant)) ./ (2*a);

        % Output
        if equation == 1
            Kd490 = Kd490_1a;
            s_Kd490 = s_Kd490_1a;
        elseif equation == 2
            Kd490 = Kd490_2a;
            s_Kd490 = s_Kd490_2a;
        end
    end
end