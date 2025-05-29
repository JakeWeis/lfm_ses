function saveOutput(root,platformID,ProfileInfo,Data)

%% CMD message: start
if isempty(getCurrentTask)
    fprintf('Saving <strong>output</strong>...');
end

%% Write everything to .mat file
s = struct('ProfileInfo',ProfileInfo,'Data',Data);
save(fullfile(root.output,[platformID(1:end-3),'_PROCESSED.mat']), '-fromstruct', s)

%% Write DATA to NetCDF file
% Open source file
ncSource = netcdf.open(fullfile(root.input,platformID), 'NOWRITE');

try
    % Get source file information
    [numDims, numVars, numGlobalAtts, ~] = netcdf.inq(ncSource);
    
    % Create new file
    if ~isfolder(fullfile(root.output,'NC'))
        mkdir(fullfile(root.output,'NC'))
    end
    ncTarget = netcdf.create(fullfile(root.output,'NC',[platformID(1:end-3),'_PROCESSED.nc']), 'NETCDF4');
    
    % Write/copy global attributes
    netcdf.putAtt(ncTarget, netcdf.getConstant('NC_GLOBAL'), 'POSTPROCESS', sprintf('J. WEIS, DATE: %s, https://github.com/JakeWeis/lfm_ses',char(datetime('now'))));
    netcdf.putAtt(ncTarget, netcdf.getConstant('NC_GLOBAL'), 'ORIGINAL_FILE_NAME', sprintf('%s (global attributes below copied from original file)',platformID));
    for iAtt = 0 : numGlobalAtts-1
        attname = netcdf.inqAttName(ncSource, netcdf.getConstant('NC_GLOBAL'), iAtt);
        attvalue = netcdf.getAtt(ncSource, netcdf.getConstant('NC_GLOBAL'), attname);
        netcdf.putAtt(ncTarget, netcdf.getConstant('NC_GLOBAL'), attname, attvalue);
    end
    
    % Define dimensions
    dimIds = containers.Map('KeyType', 'char', 'ValueType', 'double');
    dimIds('N_PROF') = netcdf.defDim(ncTarget, 'N_PROF', Data.MetaData.nProfs);
    dimIds('N_LEVELS') = netcdf.defDim(ncTarget, 'N_LEVELS', size(Data.Processed.DEPTH,1));
    
    %% Define variables
    varIds = containers.Map('KeyType', 'char', 'ValueType', 'double');
    
    % JULIAN DATE
    varname = 'JULD';
    varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', dimIds('N_PROF'));
    netcdf.putAtt(ncTarget, varIds(varname), 'long_name', 'Julian day (UTC) of the station relative to REFERENCE_DATE_TIME');
    netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', 'time');
    netcdf.putAtt(ncTarget, varIds(varname), 'units', 'days since 1950-01-01 00:00:00 UTC');
    netcdf.putAtt(ncTarget, varIds(varname), 'conventions', 'Relative julian days with decimal part (as parts of day)');
    netcdf.putAtt(ncTarget, varIds(varname), 'axis', 'T');
    netcdf.putAtt(ncTarget, varIds(varname), 'resolution', '1e-05');
    
    % LATITUDE
    varname = 'LATITUDE';
    varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', dimIds('N_PROF'));
    netcdf.putAtt(ncTarget, varIds(varname), 'long_name', 'Latitude of the station, best estimate');
    netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', varname);
    netcdf.putAtt(ncTarget, varIds(varname), 'units', 'degree_north');
    netcdf.putAtt(ncTarget, varIds(varname), 'valid_min', -90);
    netcdf.putAtt(ncTarget, varIds(varname), 'valid_max', 90);
    netcdf.putAtt(ncTarget, varIds(varname), 'axis', 'Y');

    % LONGITUDE
    varname = 'LONGITUDE';
    varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', dimIds('N_PROF'));
    netcdf.putAtt(ncTarget, varIds(varname), 'long_name', 'Longitude of the station, best estimate');
    netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', varname);
    netcdf.putAtt(ncTarget, varIds(varname), 'units', 'degree_east');
    netcdf.putAtt(ncTarget, varIds(varname), 'valid_min', -90);
    netcdf.putAtt(ncTarget, varIds(varname), 'valid_max', 90);
    netcdf.putAtt(ncTarget, varIds(varname), 'axis', 'Y');

    % PRESSURE
    varname = 'PRES';
    varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
    netcdf.putAtt(ncTarget, varIds(varname), 'long_name', 'Sea water pressure, equals 0 at sea-level');
    netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', 'sea_water_pressure');
    netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'regularised (1-m depth grid)');
    netcdf.putAtt(ncTarget, varIds(varname), 'units', 'decibar');
    netcdf.putAtt(ncTarget, varIds(varname), 'valid_min', 0);
    netcdf.putAtt(ncTarget, varIds(varname), 'valid_max', 12000);
    netcdf.putAtt(ncTarget, varIds(varname), 'resolution', 0.001);
    netcdf.putAtt(ncTarget, varIds(varname), 'axis', 'Z');

    % DEPTH
    varname = 'DEPTH';
    varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
    netcdf.putAtt(ncTarget, varIds(varname), 'long_name', 'Depth below sea level (derived from pressure and latitude using GSW equations)');
    netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', 'depth_below_sl');
    netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'regularised (1-m depth grid)');
    netcdf.putAtt(ncTarget, varIds(varname), 'units', 'm');
    netcdf.putAtt(ncTarget, varIds(varname), 'valid_min', -12000);
    netcdf.putAtt(ncTarget, varIds(varname), 'valid_max', 0);
    netcdf.putAtt(ncTarget, varIds(varname), 'resolution', 0.001);
    netcdf.putAtt(ncTarget, varIds(varname), 'axis', 'Z');

    % TEMPERATURE
    varname = 'TEMP';
    varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
    netcdf.putAtt(ncTarget, varIds(varname), 'long_name', 'Sea water temperature in-situ ITS-90 scale');
    netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', 'temperature');
    netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'regularised (1-m depth grid)');
    netcdf.putAtt(ncTarget, varIds(varname), 'units', 'degree_Celsius');
    netcdf.putAtt(ncTarget, varIds(varname), 'valid_min', -2.5);
    netcdf.putAtt(ncTarget, varIds(varname), 'valid_max', 40);
    netcdf.putAtt(ncTarget, varIds(varname), 'resolution', 0.001);

    % SALINITY
    varname = 'PSAL';
    varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
    netcdf.putAtt(ncTarget, varIds(varname), 'long_name', 'Sea water practical salinity');
    netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', 'practical_salinity');
    netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'regularised (1-m depth grid)');
    netcdf.putAtt(ncTarget, varIds(varname), 'units', 'psu');
    netcdf.putAtt(ncTarget, varIds(varname), 'valid_min', 2);
    netcdf.putAtt(ncTarget, varIds(varname), 'valid_max', 41);
    netcdf.putAtt(ncTarget, varIds(varname), 'resolution', 0.0001);

    % FLUORESCENCE
    if isfield(ProfileInfo,'FLUO')
        if ~all(ProfileInfo.FLUO.noData)
            % Regularised
            varname = 'FLUO';
            varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
            netcdf.putAtt(ncTarget, varIds(varname), 'long_name', 'Chlorophyll fluorescence');
            netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', 'mass_concentration_of_chlorophyll_a_in_sea_water');
            netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'regularised (1-m depth grid)');
            netcdf.putAtt(ncTarget, varIds(varname), 'units', 'mg/m3');

            % Dark-corrected
            varname = 'FLUO_DRK';
            varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
            netcdf.putAtt(ncTarget, varIds(varname), 'long_name', 'Chlorophyll fluorescence');
            netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', 'mass_concentration_of_chlorophyll_a_in_sea_water');
            netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'regularised, dark-corrected');
            netcdf.putAtt(ncTarget, varIds(varname), 'units', 'mg/m3');

            % NPQ-corrected
            varname = 'FLUO_DRK_NPQ';
            varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
            netcdf.putAtt(ncTarget, varIds(varname), 'long_name', 'Chlorophyll fluorescence');
            netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', 'mass_concentration_of_chlorophyll_a_in_sea_water');
            netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'regularised, dark-corrected, NPQ-corrected');
            netcdf.putAtt(ncTarget, varIds(varname), 'units', 'mg/m3');

            % Smoothed
            varname = 'FLUO_DRK_NPQ_FIT';
            varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
            netcdf.putAtt(ncTarget, varIds(varname), 'long_name', 'Chlorophyll fluorescence');
            netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', 'mass_concentration_of_chlorophyll_a_in_sea_water');
            netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'regularised, dark-corrected, NPQ-corrected, smoothed (bspline fit)');
            netcdf.putAtt(ncTarget, varIds(varname), 'units', 'mg/m3');
        end
    end

    % LIGHT
    struct_varnames = {'LIGHT','PAR','IRR490'};
    while ~isempty(struct_varnames)
        if isfield(Data.Processed,struct_varnames{1})
            % Specify attributes depending on light variable
            switch struct_varnames{1}
                case 'LIGHT'
                    varname_0 = 'LIGHT';
                    att_long_name = 'Photosynthetic photon flux density (PPFD)';
                    att_standard_name = 'PPFD';
                    att_units = 'µmol/m2/s';
                case 'PAR'
                    varname_0 = 'DOWNWELLING_PAR';
                    att_long_name = 'Downwelling photosynthetically available radiation';
                    att_standard_name = 'downwelling_PAR_in_sea_water';
                    att_units = 'µmol/m2/s';
                case 'IRR490'
                    varname_0 = 'DOWN_IRRADIANCE490';
                    att_long_name = 'Downwelling irradiance at 490 nanometers';
                    att_standard_name = 'downwelling_irradiance_490nm';
                    att_units = 'µmol/m2/s';
            end

            % Define variables
            % Regularised
            varname = varname_0;
            varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
            netcdf.putAtt(ncTarget, varIds(varname), 'long_name', att_long_name);
            netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', att_standard_name);
            netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'regularised (1-m depth grid)');
            netcdf.putAtt(ncTarget, varIds(varname), 'units', att_units);

            % Dark-corrected
            varname = [varname_0,'_DRK'];
            varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
            netcdf.putAtt(ncTarget, varIds(varname), 'long_name', att_long_name);
            netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', att_standard_name);
            netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'regularised, dark-corrected');
            netcdf.putAtt(ncTarget, varIds(varname), 'units', att_units);

            % Saturation-corrected
            varname = [varname_0,'_DRK_SAT'];
            varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
            netcdf.putAtt(ncTarget, varIds(varname), 'long_name', att_long_name);
            netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', att_standard_name);
            netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'regularised, dark-corrected, saturation-corrected');
            netcdf.putAtt(ncTarget, varIds(varname), 'units', att_units);

            % Smoothed
            varname = [varname_0,'_DRK_SAT_FIT'];
            varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
            netcdf.putAtt(ncTarget, varIds(varname), 'long_name', att_long_name);
            netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', att_standard_name);
            netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'regularised, dark-corrected, saturation-corrected, smoothed (bspline fit)');
            netcdf.putAtt(ncTarget, varIds(varname), 'units', att_units);

            % Surface extrapolated
            varname = [varname_0,'_DRK_SAT_FIT_SURF'];
            varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
            netcdf.putAtt(ncTarget, varIds(varname), 'long_name', att_long_name);
            netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', att_standard_name);
            netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'regularised, dark-corrected, saturation-corrected, smoothed (bspline fit), extrapolated to surface');
            netcdf.putAtt(ncTarget, varIds(varname), 'units', att_units);

            % Light attenuation coefficient (derived from
            varname = ['KD_',varname_0];
            varIds(varname) = netcdf.defVar(ncTarget, varname, 'double', [dimIds('N_LEVELS'),dimIds('N_PROF')]);
            netcdf.putAtt(ncTarget, varIds(varname), 'long_name', 'Light attenuation coefficient');
            netcdf.putAtt(ncTarget, varIds(varname), 'standard_name', 'light_attenuation_coefficient');
            netcdf.putAtt(ncTarget, varIds(varname), 'POST_PROCESSING', 'attenuation derived from PPFD');
            netcdf.putAtt(ncTarget, varIds(varname), 'units', '1/m');
            
            % Remove varname from cell variable
            struct_varnames(1) = [];
        else
            % Remove varname from cell variable
            struct_varnames(1) = [];
        end
    end
    
    % End define mode
    netcdf.endDef(ncTarget);

    %% Write data to variables
    netcdf.putVar(ncTarget, varIds('JULD'), convertTo(ProfileInfo.General.Date,'juliandate'));
    netcdf.putVar(ncTarget, varIds('LATITUDE'), ProfileInfo.General.Lat);
    netcdf.putVar(ncTarget, varIds('LONGITUDE'), ProfileInfo.General.Lon);
    netcdf.putVar(ncTarget, varIds('PRES'), Data.Processed.PRES.Reg);
    netcdf.putVar(ncTarget, varIds('DEPTH'), Data.Processed.DEPTH);
    netcdf.putVar(ncTarget, varIds('TEMP'), Data.Processed.TEMP.Reg);
    netcdf.putVar(ncTarget, varIds('PSAL'), Data.Processed.PSAL.Reg);
    if isKey(varIds,'FLUO')
        netcdf.putVar(ncTarget, varIds('FLUO'), Data.Processed.FLUO.Reg);
        netcdf.putVar(ncTarget, varIds('FLUO_DRK'), Data.Processed.FLUO.RegDrk);
        netcdf.putVar(ncTarget, varIds('FLUO_DRK_NPQ'), Data.Processed.FLUO.RegDrkNPQ);
        netcdf.putVar(ncTarget, varIds('FLUO_DRK_NPQ_FIT'), Data.Processed.FLUO.RegDrkNPQFitAll);
    end
    if isKey(varIds,'LIGHT')
        netcdf.putVar(ncTarget, varIds('LIGHT'), Data.Processed.LIGHT.lin.Reg);
        netcdf.putVar(ncTarget, varIds('LIGHT_DRK'), Data.Processed.LIGHT.lin.RegDrk);
        netcdf.putVar(ncTarget, varIds('LIGHT_DRK_SAT'), Data.Processed.LIGHT.lin.RegDrkSat);
        netcdf.putVar(ncTarget, varIds('LIGHT_DRK_SAT_FIT'), Data.Processed.LIGHT.lin.RegDrkSatFitAll);
        netcdf.putVar(ncTarget, varIds('LIGHT_DRK_SAT_FIT_SURF'), Data.Processed.LIGHT.lin.RegDrkSatFitSurf);
        netcdf.putVar(ncTarget, varIds('KD_LIGHT'), Data.Processed.LIGHT.Kd.FitAll);
    end
    if isKey(varIds,'DOWNWELLING_PAR')
        netcdf.putVar(ncTarget, varIds('DOWNWELLING_PAR'), Data.Processed.DOWNWELLING_PAR.lin.Reg);
        netcdf.putVar(ncTarget, varIds('DOWNWELLING_PAR_DRK'), Data.Processed.DOWNWELLING_PAR.lin.RegDrk);
        netcdf.putVar(ncTarget, varIds('DOWNWELLING_PAR_DRK_SAT'), Data.Processed.DOWNWELLING_PAR.lin.RegDrkSat);
        netcdf.putVar(ncTarget, varIds('DOWNWELLING_PAR_DRK_SAT_FIT'), Data.Processed.DOWNWELLING_PAR.lin.RegDrkSatFitAll);
        netcdf.putVar(ncTarget, varIds('DOWNWELLING_PAR_DRK_SAT_FIT_SURF'), Data.Processed.DOWNWELLING_PAR.lin.RegDrkSatFitSURF);
        netcdf.putVar(ncTarget, varIds('KD_PAR'), Data.Processed.DOWNWELLING_PAR.Kd.FitAll);
    end
    if isKey(varIds,'DOWN_IRRADIANCE490')
        netcdf.putVar(ncTarget, varIds('IRR490'), Data.Processed.DOWN_IRRADIANCE490.lin.Reg);
        netcdf.putVar(ncTarget, varIds('IRR490_DRK'), Data.Processed.DOWN_IRRADIANCE490.lin.RegDrk);
        netcdf.putVar(ncTarget, varIds('IRR490_DRK_SAT'), Data.Processed.DOWN_IRRADIANCE490.lin.RegDrkSat);
        netcdf.putVar(ncTarget, varIds('IRR490_DRK_SAT_FIT'), Data.Processed.DOWN_IRRADIANCE490.lin.RegDrkSatFitAll);
        netcdf.putVar(ncTarget, varIds('IRR490_DRK_SAT_FIT_SURF'), Data.Processed.DOWN_IRRADIANCE490.lin.RegDrkSatFitSurf);
        netcdf.putVar(ncTarget, varIds('KD_490'), Data.Processed.DOWN_IRRADIANCE490.Kd.FitAll);
    end

    % Close files
    netcdf.close(ncSource);
    netcdf.close(ncTarget);    
    
catch ME
    % Clean up in case of error
    if exist('ncSource', 'var')
        netcdf.close(ncSource);
    end
    if exist('ncTarget', 'var')
        netcdf.close(ncTarget);
    end
    rethrow(ME);
end

%% Write ProfileInfo to excel sheet
table_names = fieldnames(ProfileInfo);
for iT = 1 : numel(table_names)
    writetable(ProfileInfo.(table_names{iT}),fullfile(root.output,'NC',[platformID(1:end-3),'_ProfileInfo.xlsx']),'Sheet',table_names{iT})
end

%% CMD message: done
if isempty(getCurrentTask)
    fprintf('\b\b \x2713\n')
end

end
