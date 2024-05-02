clc

%% Project root/input/output directories
% Get root directory and add to MATLAB path
root.proj = mfilename('fullpath');
i_filesep = strfind(root.proj,filesep);
root.proj(i_filesep(end-1)+1:end) = [];
addpath(genpath(root.proj))

root.data.seal      = [root.proj '00_data' filesep '01_seal' filesep];
root.data.float     = [root.proj '00_data' filesep '02_float' filesep];

% input data directory (TO BE SPECIFIED AS INPUT TO)
% root.input = '/Volumes/PhData/PD DATA/MEOP-CTD_2024-03-08/SUBSET/FLUO_LIGHT';
root.input          = '/Volumes/PhData/PD DATA/MEOP-CTD_2024-03-08/SUBSET/FLUO_LIGHT';
if ~strcmp(root.input(end),filesep)
    root.input      = [root.input filesep];
end
root.output         = [root.input 'processed_output_NEW' filesep];
if ~isfolder(root.output)
    mkdir(root.output)
end

%% Setting default processing parameters and loading ETOPO 2022 bathymetry data
defaultPars = setDefaults(root);

if ~exist('bathymetry','var')
    fprintf('Loading <strong>ETOPO 2022 Global Relief Model</strong>...');
    [bathymetry.data, bathymetry.ref] = readgeoraster([root.proj, '00_data', filesep, '00_etopo', filesep, 'ETOPO_2022_v1_60s_N90W180_bed.tif']);
    bathymetry.data = flipud(bathymetry.data);
    bathymetry.lon = bathymetry.ref.LongitudeLimits(1) : bathymetry.ref.CellExtentInLongitude : bathymetry.ref.LongitudeLimits(2) - bathymetry.ref.CellExtentInLongitude;
    bathymetry.lat = bathymetry.ref.LatitudeLimits(1) : bathymetry.ref.CellExtentInLatitude : bathymetry.ref.LatitudeLimits(2) - bathymetry.ref.CellExtentInLatitude;    
    fprintf('\b\b \x2713\n')
else
    fprintf('Loading <strong>ETOPO 2022 Global Relief Model</strong>. \x2713\n');
end

%% Tag data processing
fprintf('\n<strong>Begin processing tags.</strong>\n\n')
% Iterate through all NetCDF files in the input directory
allFiles = dirPaths([root.input '*.nc']);
nTags = numel(allFiles);
for iTag = 1 : numel(allFiles)
    %% Processing
    % Get/display tag name currently being processed
    tagRef = allFiles(iTag).name;
    disp(repmat('-',1,numel(['Processing tag ',num2str(iTag),'/',num2str(nTags),': ',tagRef])))
    fprintf(['Processing tag %01.0f/%01.0f: <strong>',tagRef,'</strong>\n'],iTag,nTags);
    disp(repmat('-',1,numel(['Processing tag ',num2str(iTag),'/',num2str(nTags),': ',tagRef])))
    
    % Load data
    [tagData,tagMetadata,tagProcessed,ProfileInfo] = loadData(root,tagRef,bathymetry,defaultPars);
    

    % Process PAR data
    [tagProcessed,ProfileInfo_PAR] = processPAR(tagMetadata,tagProcessed,ProfileInfo,defaultPars);

    % Process Fluorescence data
    [tagProcessed,ProfileInfo_FLUO] = processFLUO(tagMetadata,tagProcessed,ProfileInfo,ProfileInfo_PAR,defaultPars);

    %% Save output
    save([root.output,tagRef(1:end-3),'_PROCESSED.mat'], ...
        'ProfileInfo','ProfileInfo_PAR','ProfileInfo_FLUO', ...
        'tagData','tagMetadata','tagProcessed')
       
end