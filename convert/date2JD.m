function JD = date2JD(year, month, day, hour, minute, second)
% JD = DATE2JD(year, month, day, hour, minute, second) converts the 
% given calendar date to the corresponding Julian date. Leapseconds 
% are not taken into account, but leapyears are. 
%
% All arguments may be vectors or matrices, as long as they all 
% have the same amount of elements. In case the arguments are non-
% scalar, the output will have the same size as the year argument.

% Authors
% .�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 15/Jun/2009 (Rody)

% (calculation simply copied from http://en.wikipedia.org/wiki/Julian_date)

    % parse input
    zero = zeros(size(year));
    if (nargin == 3), [hour, minute, second] = deal(zero); end
    if isempty(hour),   hour   = zero; end
    if isempty(minute), minute = zero; end
    if isempty(second), second = zero; end

    % leapyear corrections
    JanFeb        = month <= 2;
    year(JanFeb)  = year(JanFeb) - 1;
    month(JanFeb) = month(JanFeb) + 12;

    % Julian date
    JD = floor(365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
         floor( year/100.0 ) + floor(floor(year/100.0)/4.0) + ...   % leap year corrections
         day - 1524.5 + (hour + minute/60 + second/3600)/24;        % days and fraction of day                                 
end
