function defaultPars = setDefaults(root)
% Setting default values

fprintf('Setting default processing parameters...');

%% default parameters (MISC)
% Depth vector (m) onto which data will be interpolated, ensure that depths are a 1-by-n array of negative values
defaultPars.depthInterpGrid = -(1 : 1000);
defaultPars.depthInterpGrid = -abs(reshape(defaultPars.depthInterpGrid,numel(defaultPars.depthInterpGrid),1));

% minimum upper/lower boundary required for profiles to be imported in linear functional model (must be negative)
defaultPars.upperDepthBound = -5;
defaultPars.lowerDepthBound = -200;
defaultPars.upperDepthBound = -abs(defaultPars.upperDepthBound);
defaultPars.lowerDepthBound = -abs(defaultPars.lowerDepthBound);

% minimum profile depth (m), must be negative
defaultPars.minProfileDepth = -20;
defaultPars.minProfileDepth = -abs(defaultPars.minProfileDepth); % Ensure that depth is negative
defaultPars.minProfileLength = 20;          % minimum profile length (m between first and last obs)

% set min bathy threshold for Chl-a/light model to be run on profile
% (case 1 waters e.g. open ocean)
defaultPars.bathyMinModelThreshold = -1500; % (meters)
% criteria relative to profile depth to consider open ocean vs coastal
% when no bathymetry data available (ex. when no loc)
defaultPars.minProfDepthOpenOcean = - 500; % (meters)

defaultPars.doPlots = 'YES';                % DO/DO NOT display plots in preprocessing
defaultPars.tsDiagDepth = 200;              % depth at which pick up T and S values for T S Diagram

%% default parameters (CHLA)
defaultPars.CHL.maxCHLA_depth = -200;       % define depth at which to consider dark signal (CHLA "absolute zero")
defaultPars.CHL.maxCHLA_depth = -abs(defaultPars.CHL.maxCHLA_depth); % Ensure that depth is negative
defaultPars.CHL.delta_for_cstCHLA = 0.01;   % define standard delta between min and max of profile to eliminate profile for being a constant profile

%% default parameters (PAR)
% use of solar position function
defaultPars.PAR.UTC_offset = 0;             % [hrs] offset from UTC, during standard time
% Date and time values are always in universal time coordinates (UTC)
% see seamammal_user_manual_version1.0.pdf
defaultPars.PAR.DST = false;                % local time is without daylight savings time

% daylight filter
defaultPars.PAR.SolarAlt_DayTime = 0;	    % threshold to be applied for daylight filter (>0º solar elevation)
defaultPars.PAR.delta_constantPAR = 1E-8;   % define standard delta between min and max of profile to eliminate profile for being a constant profile
defaultPars.PAR.Epsilon = 1.0000e-6;	    % set Epsilon>0 for PAR negative values (further use of log)

% dark correction
% select n deepest per tag for dark correction
defaultPars.PAR.nDarkProfiles = 10;

%% default parameters (functional data analysis)
defaultPars.LFM = load(fullfile(root.data.seal, 'LFMcore.mat'),'LFMcore');
defaultPars.LFM.nBasis = 30;                  % number of basis functions in functional space
defaultPars.LFM.nOrder = 4;                   % order of B-splines
defaultPars.LFM.topBound_lfm = 5;
defaultPars.LFM.botBound_lfm = 200;
defaultPars.LFM.lambda = 0.03;              % smoothing parameter for functional fit
defaultPars.LFM.mybL = defaultPars.LFM.LFMcore.mybL;
defaultPars.LFM.pp = (defaultPars.LFM.topBound_lfm : defaultPars.LFM.botBound_lfm).';

%% CMD message: done
fprintf('\b\b \x2713\n')

end