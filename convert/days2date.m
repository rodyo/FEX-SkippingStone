function [year, month, day, hours, minute, second] = days2date(days)
% DAYS2DATE         
%
% Convert an amount of days since 1/1/2000 12:00 to calendar dates. 
    
    % Julian date at 2000/1/1 12:00
    JD = 2.4515445e6;    
    % convert to date
    [year, month, day, hours, minute, second] = JD2date(JD + days);
    
end
