function defVals = setDefaults(root)
% Setting default values

fprintf('Setting default parameters...');

%% default parameters (MISC)
% Depth vector (m) onto which data will be interpolated, ensure that depths are a 1-by-n array of negative values
defVals.depthInterpGrid = -(1 : 1000);
defVals.depthInterpGrid = -abs(reshape(defVals.depthInterpGrid,numel(defVals.depthInterpGrid),1));

% minimum upper/lower boundary required for profiles to be imported in linear functional model (must be negative)
defVals.upperDepthBound = -50;
defVals.lowerDepthBound = -200;
defVals.upperDepthBound = -abs(defVals.upperDepthBound);
defVals.lowerDepthBound = -abs(defVals.lowerDepthBound);

% minimum profile depth (m), must be negative
defVals.minProfileDepth = -20;
defVals.minProfileDepth = -abs(defVals.minProfileDepth); % Ensure that depth is negative
defVals.minProfileLength = 20;          % minimum profile length (m between first and last obs)

% set min bathy threshold for Chl-a/light model to be run on profile
% (case 1 waters e.g. open ocean)
defVals.bathyMinModelThreshold = -1500; % (meters)
% criteria relative to profile depth to consider open ocean vs coastal
% when no bathymetry data available (ex. when no loc)
defVals.minProfDepthOpenOcean = - 500; % (meters)

defVals.doPlots = 'YES';                % DO/DO NOT display plots in preprocessing
defVals.tsDiagDepth = 200;              % depth at which pick up T and S values for T S Diagram

%% default parameters (CHLA)
defVals.CHL.maxCHLA_depth = -200;       % define depth at which to consider dark signal (CHLA "absolute zero")
defVals.CHL.maxCHLA_depth = -abs(defVals.CHL.maxCHLA_depth); % Ensure that depth is negative
defVals.CHL.delta_for_cstCHLA = 0.01;   % define standard delta between min and max of profile to eliminate profile for being a constant profile

%% default parameters (PAR)
% use of solar position function
defVals.PAR.UTC_offset = 0;             % [hrs] offset from UTC, during standard time
% Date and time values are always in universal time coordinates (UTC)
% see seamammal_user_manual_version1.0.pdf
defVals.PAR.DST = false;                % local time is without daylight savings time

% daylight filter
defVals.PAR.SolarAlt_DayTime = 20;	    % threshold to be applied for daylight filter (>20º solar elevation)
defVals.PAR.delta_constantPAR = 1E-8;   % define standard delta between min and max of profile to eliminate profile for being a constant profile
defVals.PAR.Epsilon = 1.0000e-6;	    % set Epsilon>0 for PAR negative values (further use of log)

% dark correction
% select a few profiles per day for dark correction: the nProfiles deepest per day
defVals.PAR.nProfilesPerDay_dark = 5;

%% default parameters (functional data analysis)
defVals.LFM = load([root.data.seal 'LFMcore.mat'],'LFMcore');
defVals.LFM.nBasis = 30;                  % number of basis functions in functional space
defVals.LFM.nOrder = 4;                   % order of B-splines
defVals.LFM.topBound_lfm = 5;
defVals.LFM.botBound_lfm = 200;
defVals.LFM.lambda = 0.03;              % smoothing parameter for functional fit
defVals.LFM.mybL = defVals.LFM.LFMcore.mybL;
defVals.LFM.pp = (defVals.LFM.topBound_lfm : defVals.LFM.botBound_lfm).';

%% CMD message: done
fprintf('\b\b \x2713\n')

end