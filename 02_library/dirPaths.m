function [dirout,systemfile] = dirPaths(folder,varargin)
% DIRPATHS outputs the same structure as dir, including an addtional full path field (combined folder and name)
% Second output indicates hidden files in the list (any file preceeded by a dot)

dirout = dir(folder);
systemfile = false(numel(dirout),1);

% Input parsing
p = inputParser;
addParameter(p,'SystemFiles',false)
parse(p,varargin{:})


for iF = 1 : numel(dirout)
   dirout(iF).path = [dirout(iF).folder,'/',dirout(iF).name];
   
   if ~isempty(regexp(dirout(iF).name,'^[.]','once'))
       systemfile(iF) = true;
   else
       systemfile(iF) = false;
   end
end

if p.Results.SystemFiles == 0
    dirout(systemfile) = [];
end

end