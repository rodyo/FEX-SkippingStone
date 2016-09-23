function MPs = user_MP_model(MPs)
% USR_MP_MODEL            Mission-specific modifications or additions to
%                         the minor planets dataset (MPCORB)
%
% This function allows you to modify the Minor Planets data-set (MPCORB)
% according to your mission's specifications. When this function is called,
% you will get a variable named [MPs], in which the following fields are
% defined:
%
%   MPs.as                      semi-major axes                          [km]
%   MPs.es                      eccentricities                           [-]
%   MPs.is                      inclinations                             [rad]
%   MPs.Omegas                  longitude of ascending node              [rad]
%   MPs.omegas                  argument of perihelion                   [rad]
%   MPs.E0s                     Eccentric anomalies at J2000.0           [rad]
%   MPs.M0s                     Mean anomalies at J2000.0                [rad]
%   MPs.x0s                     Cartesian coordinates at J2000.0         [km and km/s]
%   MPs.theta0s                 true anomalies at J2000.0                [rad]
%   MPs.bs                      semi-minor axes                          [km]
%   MPs.Cs                      location of the elliptic center          [km]
%   MPs.Us                      unit vector towards perihelion           [1]
%   MPs.Vs                      unit vector towards positive [b]         [1]
%   MPs.Ws                      unit normal vector to Us and Vs          [1]
%   MPs.epoch_MJD               epoch in modified Julian dates           [days]
%   MPs.epoch_days_past_J2000   epoch in days past J2000.0               [days]
%   MPs.number_of_observations  number of observations made of the MP    [#]
%   MPs.number_of_oppositions   number of oppositions of the MP          [#]
%   MPs.first_observed          year the MP was first observed           [year]
%   MPs.last_observed           year the MP was last observed            [year]
%   MPs.types                   type of minor planet (see MPCORB format) [-]
%   MPs.names                   name or readible designation number      [-]
%   MPs.numbers                 minor planet number                      [-]
%   MPs.ns                      mean motions                             [rad/s]
%   MPs.Ts                      orbital periods                          [s]
%   MPs.Hs                      angular momentum vector (r x V)          [km2/s]
%   MPs.Hnorms                  norms of angular momentum vectors        [1]
%   MPs.peris                   perihelion radii                         [km]
%   MPs.apos                    aphelion radii                           [km]
%   MPs.number_of_MPs           total number of MPs                      [#]
%
% Note that each of these fields has a number of elements equal to the
% total amount of minor planets in the data set, i.e., 
%
%   numel(MPs.as)    = 399959 = MPs.number_of_MPs,
%   numel(MPs.types) = 399959 = MPs.number_of_MPs, 
%   ...
%
% This function allows you to add or modify fields to the dataset, in case
% your mission requires information aside from the information already
% contained in MPs. 
%
% It is generally *NOT* a good idea to *REMOVE* fields; this will usually 
% result in cryptic errors during the optimization.
%
% When the information in MPs is enough for your mission, just leave this
% function empty.    
    
    
    
    % The following is mission-specific data for 
    % Rody Oldenhuis: "Trajectory optimization to the Solar bow shock and
    % several minor planets" (MSc thesis). 

    % I just need scientific values in my cost-functions 
    % (Q_{sci} in my thesis):
    
    % initialize
    MPs.scivalue = ones(size([MPs.as; MPs.number_of_observations]));    
    % all 'normal' objects get score of 1
    
    % All other Q_{sci} depend on the type
    % (see MPCORB-header for more info)
    MPs.scivalue(MPs.types == 2 ) = 2;          % Aten
    MPs.scivalue(MPs.types == 3 ) = 2;          % Apollo
    MPs.scivalue(MPs.types == 5 ) = 2;          % Object with q < 1.381 AU
    MPs.scivalue(MPs.types == 6 ) = 2;          % Object with q < 1.523 AU
    MPs.scivalue(MPs.types == 7 ) = 2;          % Object with q < 1.665 AU
    MPs.scivalue(MPs.types == 8 ) = 2;          % Hilda
    MPs.scivalue(MPs.types == 9 ) = 8;          % Jupiter Trojan
    MPs.scivalue(MPs.types == 10) = 10;         % Centaur
    MPs.scivalue(MPs.types == 14) = 10;         % Plutino
    MPs.scivalue(MPs.types == 15) = 10;         % Other resonant TNO
    MPs.scivalue(MPs.types == 16) = 10;         % Cubewano
    MPs.scivalue(MPs.types == 17) = 20;         % Scattered disk
    MPs.scivalue(MPs.types == 32768) = 100;     % Object is PHA            
            
end % user MP model
