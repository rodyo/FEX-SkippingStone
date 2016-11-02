function [year, month, day, hour, minute, second] = MJD2date(MJD)
% MJD2DATE              Convert modified Julian dates to Gregorian
%                       calendar dates
%
% [year, month, day, hour, minute, second] = MJD2DATE(MJD) converts
% the given modified Julian date to the Gregorian calendar date.
% Works correctly for scalar/vector/matrix input.
%
% See also JD2date, days2date, date2MJD, date2JD.

% Author:
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com

% If you find this work useful, please consider a small donation:
% https://www.paypal.me/RodyO/3.5

    % convert to JD
    JD = MJD + 2400000.5;
    % convert JD to Gregorian date
    [year, month, day, hour, minute, second] = JD2date(JD);

end
