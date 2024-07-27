function [Fronts,Zones] = SOFronts(varargin)
% Extract Orsi, AH., Harris, U. (2019) Southern Ocean front coordinates. Coordinates are centered around the 0˚ meridian.
%
% [fronts,zones] = SOFronts;
% Returns standard Orsi fronts and zones:
% >> NBdy       Northern Southern Ocean boundary (30˚S)
%    STZ        Subtropical zone
% >> STF        Subtropical front
%    SAZ        Subantarctic Zone
% >> SAF        Subantarctic front
%    PFZ        Polar frontal zone
% >> PF         Polar front
%    AZ         Antarctic Zone
% >> SACCF      Southern Antarctic Circumpolar Current front
%    SZ         Southern zone
% >> SBdy       Southern boundary
%    SPR        Subpolar region
% 
% Add custom fronts (multiple fronts as n-by-3 cell array):
% 'add_fronts'          {'front_name',longitudes,latitudes}
%   front_name          Name of the custom front (NB: used as structure field name)
%   longitudes          Longitudes of the custom front
%   latitudes           Latitudes of the custom front
% 
% Add custom zone (multiple zones as n-by-3 cell array)
% 'add_zone'            {'zone_name','N_front','S_front'}
%   zone name           Name of the custom zone  (NB: used as structure field name)
%   N_front             Name of the northern front  (either one of the standard Orsi fronts or a custom front)
%   S_front             Name of the southern front  (either one of the standard Orsi fronts or a custom front)
%
% Example:
% front1 = {'lat_53S',linspace(-180,180,500),repmat(-53,500,1)};
% front2 = {'lat_63S',linspace(-180,180,500),repmat(-63,500,1)};
% custom_fronts = [front1;front2];
% zone1 = {'STF_53S','STF','lat_53S'};
% zone2 = {'SAF_63S','SAF','lat_63S'};
% custom_zones = [zone1,zone2];
%
% [fronts,zones] = SOFronts('add_fronts',custom_fronts,'add_zones',custom_zones);
%
% Reference: Orsi, AH., Harris, U. (2019) Fronts of the Antarctic Circumpolar Current GIS data, Ver. 1, Australian Antarctic
% Data Centre https://data.aad.gov.au/metadata/records/antarctic_circumpolar_current_fronts Accessed: 2020-11-19

% Parse input
p = inputParser;
addParameter(p,'add_fronts',{})
addParameter(p,'add_zones',{})
parse(p,varargin{:});

opt.AddFronts = p.Results.add_fronts;
opt.AddZones = p.Results.add_zones;

frontstable = readtable('/Users/jweis/MATLAB-Drive/misc/Analyses/SOFronts/OrsiFronts.csv');

% SO front lat/lon data
front_name = {'NBdy','STF','SAF','PF','SACCF','SBdy'};
for iF = 1 : length(front_name)
    % Extract from Orsi et al data table
    Fronts.(front_name{iF}).lon = frontstable{strcmp(frontstable{:,3},front_name{iF}),1};
    Fronts.(front_name{iF}).lat = frontstable{strcmp(frontstable{:,3},front_name{iF}),2};

    % Extend vectors to -180 and 180˚ lon (repeating first and last latitude value)
    Fronts.(front_name{iF}).lon = [-180;Fronts.(front_name{iF}).lon;180];
    Fronts.(front_name{iF}).lat = [Fronts.(front_name{iF}).lat(1);Fronts.(front_name{iF}).lat;Fronts.(front_name{iF}).lat(end)];
end

% SO frontal zone polygons
zone_name = {'STZ','SAZ','PFZ','AZ','SZ','SPR'};
for iZ = 1 : numel(zone_name)
    if iZ < numel(zone_name)
        Zones.(zone_name{iZ}).lon = [Fronts.(front_name{iZ}).lon;flipud(Fronts.(front_name{iZ+1}).lon)];
        Zones.(zone_name{iZ}).lat = [Fronts.(front_name{iZ}).lat;flipud(Fronts.(front_name{iZ+1}).lat)];
    else
        Zones.(zone_name{iZ}).lon = [Fronts.(front_name{iZ}).lon;linspace(180,-180,100)'];
        Zones.(zone_name{iZ}).lat = [Fronts.(front_name{iZ}).lat;repmat(-89,100,1)];
    end
end

% Add custom fronts
if ~isempty(opt.AddFronts)
    for iCF = 1 : size(opt.AddFronts,1)
        Fronts.(opt.AddFronts{iCF,1}).lon = reshape(opt.AddFronts{iCF,2},numel(opt.AddFronts{iCF,2}),1);
        Fronts.(opt.AddFronts{iCF,1}).lat = reshape(opt.AddFronts{iCF,3},numel(opt.AddFronts{iCF,3}),1);
    end
end

% Add custom zones
if ~isempty(opt.AddZones)
    for iCZ = 1 : size(opt.AddZones,1)
        Zones.(opt.AddZones{iCZ,1}).lon = [Fronts.(opt.AddZones{iCZ,2}).lon;flipud(Fronts.(opt.AddZones{iCZ,3}).lon)];
        Zones.(opt.AddZones{iCZ,1}).lat = [Fronts.(opt.AddZones{iCZ,2}).lat;flipud(Fronts.(opt.AddZones{iCZ,3}).lat)];
    end
end

end