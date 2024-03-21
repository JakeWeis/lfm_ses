% https://www.mathworks.com/matlabcentral/fileexchange/15585-datepart

function varargout = datepart(din,varargin)
% DATEPART Extract decimal year/month/day/hour/min/sec from a datenum.
%   X = DATEPART(DATE,'part') returns the date 'part' (year, month,
%   day, hour, minute or second) of the datenum date with the decimal
%   equivalent of the trailing date parts in the units of 'part'. DATE can
%   be scalar, vector or array and X will be of the some size and shape.
% 
%   [X1,X2,...] = DATEPART(DATE,'part1','part2',...) returns each decimal
%   value corresponding to each part input argument.
% 
%   Acceptable string specifiers for 'part' are:
%   Year: 'year','yr','years','yrs'
%   Month: 'month','mn','months','mns'
%   Day: 'day','days'
%   Hour: 'hour','hours','hr','hrs'
%   Minute: 'minutes','minutes','min','mins'
%   Second: 'second','seconds','sec','secs'
% 
%   Examples:
%   >> A = datenum('14-Jul-2007 15:05:12')
%   >> datepart(A,'min')
%   ans =
%        5.2000
%   >> floor(datepart(A,'min'))
%   ans =
%       5
%   >> format bank
%   >> datepart(A,'year')
%   ans =
%       2007.53
%   >> [day,yr] = datepart(A,'day','year')
%   day =
%          14.63
%   year =
%        2007.53
%
% Ian M. Howat, Applied Physics Lab, University of Washington
% ihowat@apl.washington.edu
% Version 1: 14-Jul-2007 16:11:57
%   Revision 1: 18-Jul-2007 16:48:31
%       Added multiple-input string functionality.
%%
[r,c]=size(din);
din = datevec(din(:));
if nargin < 2
    error('Not enough input arguments.')
end
varargout = cell(length(varargin),1);
for k=1:length(varargin)
   
    switch lower(varargin{k})
        case {'year','yr','years','yrs'}
            n = 1;
        case {'month','mn','months','mns'}
            n = 2;
        case {'day','days'}
            n = 3;
        case {'hour','hours','hr','hrs'}
            n = 4;
        case {'minute','minutes','min','mins'}
            n = 5;
        case {'second','seconds','sec','secs'}
            n = 6;
        otherwise
            error(['>> "',varargin{k},'" not recognized']);
    end
    
    if n == 1 %years
        %nd days in year
        nd = datenum([din(:,1)+ 1, repmat([0,0,0,0,0],[size(din,1),1])]) -...
            datenum([din(:,1), repmat([0,0,0,0,0],[size(din,1),1])]);
        dout = din(:,1)+(datenum(din) - datenum([din(:,1),repmat([1 1 0 0 0],...
            [size(din(:,1),1),1])]))./nd;
    else
        dout = din(:,6);%decimal seconds
        if n < 6; 
            dout = din(:,5) + (dout./60); %decimal minutes
        end
        if n < 5
            dout = din(:,4) + (dout./60); %decimal hours
        end
        if n < 4
            dout = din(:,3) + (dout./24); %decimal days
        end
        if n < 3 %months
            %nd days in month
            nd = datenum([din(:,1:2)+[zeros(size(din(:,1))),ones(size(din(:,1)))],...
                repmat([0,0,0,0],[size(din,1),1])]) - datenum([din(:,1:2),...
                repmat([0,0,0,0],[size(din,1),1])]);
            dout = din(:,2) + (dout./nd); % decimal months
        end
    end
    dout = reshape(dout,[r,c]);
    varargout{k} = dout;
end
