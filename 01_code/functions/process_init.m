function root = process_init(root)

%% Project root and raw data input directories
root.data.etopo  = fullfile(root.proj, '00_data','etopo');
root.data.bedmap  = fullfile(root.proj, '00_data','bedmap3');
root.data.seal   = fullfile(root.proj, '00_data','seal');
root.data.seaice  = fullfile(root.proj, '00_data','seaice');

% Check for ETOPO files and convert to binary if needed
etopo_nc_files = dirPaths(fullfile(root.data.etopo,'ETOPO_2022_v1_*_N90W180_bed.nc'));
etopo_flt_files = dirPaths(fullfile(root.data.etopo,'ETOPO_2022_v1_*_N90W180_bed.nc'));
root.missing.etopo_nc = isempty(etopo_nc_files);
root.missing.etopo_flt = isempty(etopo_flt_files);
if ~root.missing.etopo_nc
    etopo_dataset = etopo_nc_files(1).name(1:end-3);
    
    if root.missing.etopo_flt
        % If binary (.flt) files are missing, create from the NetCDF
        disp('Writing ETOPO bathymetry data to binary file (reading information from NetCDF file takes waaaaaay longer).')

        % Bathymetry
        data = ncread(fullfile(root.data.etopo,[etopo_dataset,'.nc']),'z');
        data_single = single(data);
        fid = fopen(fullfile(root.data.etopo,[etopo_dataset,'.flt']), 'w');
        fwrite(fid, data_single, 'float32');
        fclose(fid);

        % Switch logical missing file indicator
        root.missing.etopo_flt = false;
    end
else
    warning('Missing ETOPO bathymetry file. Respective bathymetry will not be determined. Refer to NOTE.txt files in ./00_data/etopo for download instructions.')
end

% Check for BEDMAP files and convert to binary if needed
bedmap_nc_files = dirPaths(fullfile(root.data.bedmap,'bedmap3.nc'));
bedmap_flt_files = dirPaths(fullfile(root.data.bedmap,'bm3_bed_topography.flt'));
root.missing.bedmap_nc = isempty(bedmap_nc_files);
root.missing.bedmap_flt = isempty(bedmap_flt_files);
if ~root.missing.bedmap_nc
    if root.missing.bedmap_flt
        % If binary (.flt) files are missing, create from the NetCDF
        disp('Writing BEDMAP bathymetry data to binary file (reading information from NetCDF file takes waaaaaay longer).')

        % Bathymetry
        data = ncread(fullfile(root.data.bedmap,'bedmap3.nc'),'bed_topography');
        data_single = single(data);
        fid = fopen(fullfile(root.data.bedmap,'bm3_bed_topography.flt'), 'w');
        fwrite(fid, data_single, 'float32');
        fclose(fid);

        % Uncertainty
        data = ncread(fullfile(root.data.bedmap,'bedmap3.nc'),'bed_uncertainty');
        data_single = single(data);
        fid = fopen(fullfile(root.data.bedmap,'bm3_bed_uncertainty.flt'), 'w');
        fwrite(fid, data_single, 'float32');
        fclose(fid);

        % Switch logical missing file indicator
        root.missing.bedmap_flt = false;
    end
else
    warning('Missing BEDMAP bathymetry file. Respective bathymetry will not be determined. Refer to NOTE.txt files in ./00_data/bedmap3 for download instructions.')
end

root.missing.sic = isempty(dirPaths(fullfile(root.data.seaice,'*.nc')));
if root.missing.sic
    warning('Missing sea ice concentration maps. No sea ice concentration will be determined. Refer to NOTE.txt files in ./00_data/seaice for download instructions.')
end

end