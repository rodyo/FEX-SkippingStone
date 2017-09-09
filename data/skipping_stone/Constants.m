classdef ...
(Sealed) ...
Constants

    properties (Constant)
        
        % Astrodynamical constants from (From NASA/JPL Horizons)
        % and some convenient constants
        G  = 6.67300e-20;   % gravitational constant    [km+3 kg-1 s-2]
        AU = 149597870.691; % Astronomical Unit         [km]
        g0 = 9.80665/1e3;   % grav. acc. at sealevel    [km s-2]
                            % (note the conversion to km)

        % seconds, minutes, hours, days, weeks, months, years
        timeunits = {1
                     60
                     3600
                     86400
                     7*86400
                     30.471*86400
                     365.25*86400};
        
    end 
    
end
