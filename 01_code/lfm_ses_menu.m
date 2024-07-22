clc

%% Project root/input/output directories
% Get root directory and add to MATLAB path
root.proj = fileparts(fileparts(mfilename('fullpath'))); % lfm_ses_menu is in a subdirectory two levels down from the root directory
addpath(genpath(root.proj))

root.data.etopo  = fullfile(root.proj, '00_data','etopo');
root.data.seal   = fullfile(root.proj, '00_data','seal');
root.data.seaice  = fullfile(root.proj, '00_data','seaice');

% input data directory
data_base_dir       = '/Volumes/PhData/PD DATA';
root_input{1}       = '/SEALTAGS/Prydz/FLUO_LIGHT';
root_input{2}       = '/SEALTAGS/Prydz/FLUO';
root_input{3}       = '/BGCARGO/FLUO_LIGHT/Profiles';

for iInput = 3
    root.input       = fullfile(data_base_dir,root_input{iInput},'RAW');
    root.output      = fullfile(data_base_dir,root_input{iInput},'PROCESSED');
    if ~isfolder(root.output)
        mkdir(root.output)
    end

    %% Setting default processing parameters and loading ETOPO 2022 bathymetry data
    defaultPars = setDefaults(root);
    pool = gcp;

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

        % Process Fluorescence data
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