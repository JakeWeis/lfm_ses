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

%% Project root and raw data input directories
% Get LFM_SES root directory and add to MATLAB path
root.proj = fileparts(fileparts(mfilename('fullpath'))); % lfm_ses_menu is in a subdirectory two levels down from the root directory
addpath(genpath(root.proj))

root.data.etopo  = fullfile(root.proj, '00_data','etopo');
root.data.bedmap  = fullfile(root.proj, '00_data','bedmap2');
root.data.seal   = fullfile(root.proj, '00_data','seal');
root.data.seaice  = fullfile(root.proj, '00_data','seaice');

% Issue warnings if bathymetry or sea ice concentration datasets are missing
missingBathyFiles = (...
    ~exist(fullfile(root.data.etopo,'ETOPO_2022_v1_30s_N90W180_bed.nc'),"file") &&...
    ~exist(fullfile(root.data.etopo,'ETOPO_2022_v1_0s_N90W180_bed.nc'),"file")) ||...
    ~exist(fullfile(root.data.bedmap,'bedmap2_bed.flt'),"file") ||...
    ~exist(fullfile(root.data.bedmap,'bedmap2_grounded_bed_uncertainty.flt'),"file");
if missingBathyFiles
    warning('Missing ETOPO and BEDMAP bathymetry maps. No bathymetry will be determined. Refer to NOTE.txt files in ./00_data/bedmap2 and ./00_data/etopo for download instructions.')
end

missingSICFiles = isempty(dirPaths(fullfile(root.data.seaice,'*.nc')));
if missingSICFiles
    warning('Missing sea ice concentration maps. No sea ice concentration will be determined. Refer to NOTE.txt files in ./00_data/seaice for download instructions.')
end

% Raw input data directory
data_base_dir       = '/Volumes/Data/Work';
% root_input{1}       = '/SEALTAGS/Prydz/ct182/v3_20250310/';
root_input{1}       = '/SEALTAGS/Prydz/FLUO_LIGHT';
root_input{2}       = '/SEALTAGS/Prydz/FLUO';

% Set default processing parameters
defaultPars = setDefaults(root);

for iInput = 1:2
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
        [Data,ProfileInfo] = loadData(root,platformID,defaultPars);

        % Process light data
        switch Data.MetaData.platform_type
            case 'sealtag'
                if isfield(Data.Processed,'LIGHT')
                    [Data,ProfileInfo] = processRadiometry(Data,ProfileInfo,defaultPars,'LIGHT');
                end
            case 'float'
                if isfield(Data.Processed,'PAR')
                    [Data,ProfileInfo] = processRadiometry(Data,ProfileInfo,defaultPars,'DOWNWELLING_PAR');
                end
                if isfield(Data.Processed,'IRR490')
                    [Data,ProfileInfo] = processRadiometry(Data,ProfileInfo,defaultPars,'DOWN_IRRADIANCE490');
                end
        end

        % Process fluorescence data
        if isfield(Data.Processed,'FLUO')
            [Data,ProfileInfo] = processFluorometry(Data,ProfileInfo,defaultPars);
        end

        %% Save output
        saveOutput(root,platformID,ProfileInfo,Data)
   end
end