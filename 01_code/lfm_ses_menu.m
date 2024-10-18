%% lfm_ses_menu
% Main script used to execute the seal tag and BGC-Argo data processing. This script loops through NetCDF files in the
% specified seal tag/BGC-Argo data directories (multiple directories can be specified and looped through if datasets are
% stored in more than one folder) and applies a number of processing and correction steps ot the data. Processed data are
% stored in the raw data sub-directory as .mat files. This script also loads the ETOPO 2022 bathymetry dataset used to get
% the bathymetry at each profile location.
%
% The processing is divided into the following subscripts:
% 1) setDefaults.m: Set default processing parameters 
% 2) loadData.m: Loads raw data and extracts metadata from seal tag or BGC-Argo NetCDF files. Raw data is interpolated onto a
% uniform 1-m depth grid.
% 3) processPAR.m: PAR data processing
% 4) processFLUO.m: Fluorescence data processing

%% Project root and raw data input directories
% Get LFM_SES root directory and add to MATLAB path
root.proj = fileparts(fileparts(mfilename('fullpath'))); % lfm_ses_menu is in a subdirectory two levels down from the root directory
addpath(genpath(root.proj))

root.data.etopo  = fullfile(root.proj, '00_data','etopo');
root.data.seal   = fullfile(root.proj, '00_data','seal');
root.data.seaice  = fullfile(root.proj, '00_data','seaice');

% Raw input data directory
data_base_dir       = '/Volumes/Data/Work';
root_input{1}       = '/SEALTAGS/Prydz/FLUO_LIGHT';
root_input{2}       = '/SEALTAGS/Prydz/FLUO';
root_input{3}       = '/BGCARGO/FLUO_LIGHT/Profiles';

% Set default processing parameters
defaultPars = setDefaults(root);

for iInput = 1 : 2
    %% Raw data input and processing output directories
    root.input       = fullfile(data_base_dir,root_input{iInput},'RAW');
    root.output      = fullfile(data_base_dir,root_input{iInput},'PROCESSED');
    if ~isfolder(root.output)
        mkdir(root.output)
    end

    %% Load ETOPO 2022 bathymetry data
    if ~exist('bathymetry','var')
        fprintf('Loading <strong>ETOPO 2022 Global Relief Model</strong>...');
        [bathymetry.data, bathymetry.ref] = readgeoraster(fullfile(root.data.etopo, 'ETOPO_2022_v1_60s_N90W180_bed.tif'));
        bathymetry.data = flipud(bathymetry.data);
        bathymetry.lon = bathymetry.ref.LongitudeLimits(1) : bathymetry.ref.CellExtentInLongitude : bathymetry.ref.LongitudeLimits(2) - bathymetry.ref.CellExtentInLongitude;
        bathymetry.lat = bathymetry.ref.LatitudeLimits(1) : bathymetry.ref.CellExtentInLatitude : bathymetry.ref.LatitudeLimits(2) - bathymetry.ref.CellExtentInLatitude;
        fprintf('\b\b \x2713\n')
    else
        fprintf('Loading <strong>ETOPO 2022 Global Relief Model</strong>. \x2713\n');
    end

    %% Seal tag/float data processing
    fprintf('\n<strong>Begin processing platform data.</strong>\n\n')
    % Iterate through all NetCDF files in the input directory
    allFiles = dirPaths(fullfile(root.input, '*.nc'));
    nFiles = numel(allFiles);

   parfor iPlatform = 1 : numel(allFiles)
        %% Processing
        % Platform ID
        platformID = allFiles(iPlatform).name;

        % Display platform being processed
        disp(repmat('-',1,numel(sprintf('Processing platform %i/%i: %s',iPlatform,nFiles,platformID))))
        fprintf('Processing platform %i/%i: <strong>%s</strong>\n',iPlatform,nFiles,platformID);
        disp(repmat('-',1,numel(sprintf('Processing platform %i/%i: %s',iPlatform,nFiles,platformID))))

        % Load data
        [Data,ProfileInfo] = loadData(root,platformID,bathymetry,defaultPars);

        % Process light data
        [Data,ProfileInfo] = processPAR(Data,ProfileInfo,defaultPars,'PAR');
        [Data,ProfileInfo] = processPAR(Data,ProfileInfo,defaultPars,'IRR490');     % downwelling irradiance 490nm, BGC-Argo only

        % Process fluorescence data
        [Data,ProfileInfo] = processFLUO(Data,ProfileInfo,defaultPars);

        %% Save output
        s = struct('ProfileInfo',ProfileInfo,'Data',Data);
        save(fullfile(root.output,[platformID(1:end-3),'_PROCESSED.mat']), '-fromstruct', s)

    end

    %% Load data
    % for iPlatform = 1 : numel(allFiles)
    %     platformID = allFiles(iPlatform).name;
    %     load(fullfile(root.output,[platformID(1:end-3),'_PROCESSED.mat']),'Data');
    %     if ~isfield(Data.Raw,'DOWN_IRRADIANCE490')
    %         sprintf('no kd490 (%s)',ProfileInfo.General.PlatformID(1,:))
    %     end
    % end
end