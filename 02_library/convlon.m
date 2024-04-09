function NewLon = convlon(Lon,varargin)
% Converting longitudes between signed (-180 to 180˚) and unsigned (0 to 360˚) format.
%
% Input parameters
% 1) Longitude (required): Number/vector/array of longitude values to be converted.
% 2) Specifier (optional): 'signed' converts from unsigned to signed, 'unsigned' converts
%        from signed to unsigned, no input automatically determines and converts between
%        formats.

if ~isempty(varargin)
    if strcmp(varargin{1}, 'signed')
        NewLon = Lon;
        NewLon(NewLon > 180) = NewLon(NewLon > 180) - 360;
    elseif strcmp(varargin{1}, 'unsigned')
        NewLon = Lon;
        NewLon(NewLon < 0) = NewLon(NewLon < 0) + 360;
        NewLon(NewLon > 360) = NewLon(NewLon > 360) - 360;
    else
        warning('Conversion specifier not recognised (''signed'', ''unsigned'' or no specifier for automatic conversion.')
    end
else
    NewLon = Lon;
    if ~isempty(NewLon(NewLon > 180))
        NewLon(NewLon > 180) = NewLon(NewLon > 180) - 360;
        %disp('Converted from unsigned to signed')
    elseif ~isempty(NewLon(NewLon < 0))
        NewLon(NewLon < 0) = NewLon(NewLon < 0) + 360;
        %disp('Converted from signed to unsigned')
    end
end
    
end