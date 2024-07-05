function ncdata = ncload_struct(fileIn, varargin)

% ncload -- Load NetCDF variables.
%  ncload('fileIn', 'var1', 'var2', ...) loads the
%   given variables of 'fileIn' into the Matlab
%   workspace of the "caller" of this routine.  If no names
%   are given, all variables are loaded.

ncdata = struct;

if ~exist(fileIn,'file')
    sprintf('error: %s does not exist\n',fileIn)
    return
end

if isempty(varargin)
    finfo = ncinfo(fileIn);
    vars={finfo.Variables.Name};
    varargin = vars; 
end

for iVar = 1 : length(varargin)
    try
        vardata = ncread(fileIn,varargin{iVar});
        if isnumeric(vardata)
            vardata = double(vardata); 
        end

        ncdata.(varargin{iVar}) = vardata;

    catch
        sprintf('error: %s is not a variable of %s\n',varargin{iVar},fileIn)
    end
end