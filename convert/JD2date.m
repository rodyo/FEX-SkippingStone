function [year, month, day, hour, minute, second] = JD2date(JD)    
% JD2DATE          Convert Julian dates to Gregorian calendar dates
%
% [year, month, day, hour, minute, second] = JD2DATE(JD) converts 
% the given Julian date [JD] to the corresponding Gregorian calendar 
% date. Works correctly for scalar/vector/matrix input.
%
% When only [year, month, day] are requested, the returned [day] is 
% a fractional day; that is, the hours/minutes/seconds are appended to the
% integer [day] as (1/24)th, (1/24/60)th or (1/24/60/60)th of a day,
% respectively. If on the other hand the hours, minutes and seconds are
% requested explicitly, the returned [day] is an integer. For example: 
%
%   day = 31, hours = 12, minutes = 30, and hours etc. are NOT requested:
%
%       output [day] = 31.52083333...
%
%   day = 31, hours = 12, minutes = 30, and hours etc. ARE requested:
%
%       output [day]     = 31
%       output [hours]   = 12
%       output [minutes] = 30
%
% LIMITATIONS (from wikipedia:)
% 
% "For dates before 1582, the resulting date components are valid only 
% in the Gregorian proleptic calendar. This is based on the Gregorian 
% calendar but extended to cover dates before its introduction, including 
% the pre-Christian era. For dates in that era (before year 1 AD), 
% astronomical year numbering is used. This includes a year zero, which 
% immediately precedes 1 AD. Astronomical year zero is 1 BC in the proleptic 
% Gregorian calendar and, in general, proleptic Gregorian year (n BC) = 
% astronomical year (Y = 1 ? n). For astronomical year Y (Y < 1), the 
% proleptic Gregorian year is (1 - Y) BC."
%
% See also MJD2date, days2date, date2MJD, date2JD.

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 01/Nov/2009    

    % integer part
    p = fix(JD + 0.5);

    % conversion factors
    % (from http://www.astro.uu.nl/~strous/AA/en/reken/juliaansedag.html)   
    s1 = p + 68569;                     n  = floor(4*s1/146097);        
    s2 = s1 - floor((146097*n + 3)/4);  i  = floor(4000*(s2 + 1)/1461001);        
    s3 = s2 - floor(1461*i/4) + 31;     q  = floor(80*s3/2447);         
    e  = s3 - floor(2447*q/80);         s4 = floor(q/11);

    % calculate month & year        
    month = q + 2 - 12*s4 ;
    year  = 100*(n - 49) + i + s4;
    
    % output of fractional days or hours/minutes/second depends on whether 
    % hours/minute/second was actually requested
    if nargout < 4
        % add fractional days to [e]
        day = e + JD - p + 0.5;
    else
        % day = [e]
        day = e;        
        % convert fractional days to seconds
        frac = (JD + 0.5 - p)*86400;        
        % convert to hours, minutes, seconds
        hour   = floor(frac/3600);   frac = frac - 3600*hour;
        minute = floor(frac/60);
        second = frac - 60*minute;
    end

end
