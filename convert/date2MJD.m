function MJD = date2MJD(year, month, day, hour, minute, second)
% MJD = DATE2MJD(year, month, day, hour, minute, second) simply 
% converts the given date to the modified Julian date. 

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 15/Jun/2009 (Rody)

    % parse input
    zero = zeros(size(year));
    if (nargin == 3), [hour, minute, second] = deal(zero); end
    if isempty(hour),   hour   = zero; end
    if isempty(minute), minute = zero; end
    if isempty(second), second = zero; end

    % get Julian date
    JD = date2JD(year, month, day, hour, minute, second);

    % subtract the constant 
    MJD = JD - 2400000.5;

end
