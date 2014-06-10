function sec = date2sec(year, month, day, hour, minute, second)
% days = DATE2SEC(year, month, day, hour, minute, second) converts the
% given calendar date into the amount of seconds since the J2000.0 epoch.

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

    % compute amount of days since J2000.0
    days = date2days(year, month, day, hour, minute, second);

    % convert to seconds
    sec = days*86400;

end
