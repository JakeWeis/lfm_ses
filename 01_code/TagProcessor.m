%% TagProcessor
% Main script used to execute the seal tag and BGC-Argo data processing. This script loops through NetCDF files in the
% specified seal tag/BGC-Argo data directories (multiple directories can be specified and looped through if datasets are
% stored in more than one folder) and applies a number of processing and correction steps ot the data. Processed data are
% stored in the raw data sub-directory as .mat files. This script also loads the ETOPO 2022 bathymetry dataset used to get
% the bathymetry at each profile location.
%
% The processing is divided into the following subscripts:
% 1) setDefaults.m: Set default processing parameters (to change defaults, edit parameters in the script)
% 2) loadData.m: Loads raw data and extracts metadata from seal tag or BGC-Argo NetCDF files. Raw data is interpolated onto a
% uniform 1-m depth grid.
% 3) processPAR.m: PAR data processing
% 4) processFLUO.m: Fluorescence data processing

% Get LFM_SES root directory and add to MATLAB path
root.proj = fileparts(fileparts(mfilename('fullpath'))); % lfm_ses_menu is in a subdirectory two levels down from the root directory
addpath(genpath(root.proj))
root = process_init(root);

% Raw input data directory
data_base_dir       = '/Volumes/Data/Work';
% root_input{1}       = '/SEALTAGS/Prydz/FLUO_LIGHT';
% root_input{2}       = '/SEALTAGS/Prydz/FLUO';
% root_input{3}       = '/SEALTAGS/Prydz/ct182/v4_20250414';
root_input{1}       = '/SEALTAGS/Leo/Australia_filtered/DATA';
root_input{2}       = '/SEALTAGS/Leo/Australia_filtered/DATA_FULL_RES';
root_input{3}       = '/SEALTAGS/Leo/France_filtered/DATA';
root_input{4}       = '/SEALTAGS/Leo/France_filtered/DATA_FULL_RES';

% Set default processing parameters
defaultPars = setDefaults(root);

for iInput = 1:4
    %% Raw data input and processing output directories
    root.input       = fullfile(data_base_dir,root_input{iInput},'RAW');
    root.output      = fullfile(data_base_dir,root_input{iInput},'PROCESSED');
    if ~isfolder(root.output)
        mkdir(root.output)
    end

    %% Seal tag/float data processing
    fprintf('\n<strong>Begin processing platform data.</strong>\n\n')
    % Iterate through all NetCDF files in the input directory
    allFiles = dirPaths(fullfile(root.input, '*.nc'));
    % Remove "_traj" files
    allFiles(contains({allFiles.name}','_traj')) = [];
    nFiles = numel(allFiles);

   for iPlatform = 1 : numel(allFiles)
        %% Processing
        % Platform ID
        platformID = allFiles(iPlatform).name;

        % Display platform being processed
        disp(repmat('-',1,numel(sprintf('Processing platform %i/%i: %s',iPlatform,nFiles,platformID))))
        fprintf('Processing platform %i/%i: <strong>%s</strong>\n',iPlatform,nFiles,platformID);
        disp(repmat('-',1,numel(sprintf('Processing platform %i/%i: %s',iPlatform,nFiles,platformID))))

        % Load data
        [Data,ProfileInfo] = loadData(root,platformID,defaultPars);
        hasLight = isfield(Data.Processed,'LIGHT') || isfield(Data.Processed,'DOWNWELLING_PAR') || isfield(Data.Processed,'DOWN_IRRADIANCE490');
        hasFluo = isfield(Data.Processed,'FLUO');

        % Process light data
        if hasLight
            if isfield(Data.Processed,'LIGHT')
                [Data,ProfileInfo] = processRadiometry(Data,ProfileInfo,defaultPars,'LIGHT');
            end
            if isfield(Data.Processed,'PAR')
                [Data,ProfileInfo] = processRadiometry(Data,ProfileInfo,defaultPars,'DOWNWELLING_PAR');
            end
            if isfield(Data.Processed,'IRR490')
                [Data,ProfileInfo] = processRadiometry(Data,ProfileInfo,defaultPars,'DOWN_IRRADIANCE490');
            end
        end

        % Process fluorescence data
        if hasFluo
            [Data,ProfileInfo] = processFluorometry(Data,ProfileInfo,defaultPars);
        end

        % Calculate photophysiological metrics
        if hasLight && hasFluo
            [Data,ProfileInfo] = PhotoPhysiology(Data,ProfileInfo,defaultPars);
        end

        %% Save output
        saveOutput(root,platformID,ProfileInfo,Data)
   end
end