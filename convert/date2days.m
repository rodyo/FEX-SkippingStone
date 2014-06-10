function days = date2days(year, month, day, hour, minute, second)
% days = DATE2DAYS(year, month, day, hour, minute, second) converts the
% given calendar date into the amount of days since the J2000.0 epoch.

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

    % Julian date at Jan 1st, 2000, noon (12:00). 
    null = 2451545.0;         

    % Julian date of reqested epoch
    JD = date2JD(year, month, day, hour, minute, second);        

    % number of days since the reference
    days = (JD - null);                                     

end
