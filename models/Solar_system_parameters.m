function [model, constants] = Solar_system_parameters(model, constants)
% SOLAR_SYSTEM_PARAMETERS   Definitions of various parameters for all the planets.
%

% TODO: add comments (which parameters, etc)

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com
% Affiliation: Delft University of Technology
%              Partly performed at ISAS/JAXA

    % Last Edited: 29/Oct/2010
        
    % name of central body    
    model.CentralBody{1} = 'Sun';        % name
    model.CentralBody{2} = 1;            % index
    model.CentralBody{3} = 132712439940; % std. grav. parameter [km3 s-2]
    
    % convenient copy of std. grav. parameter Sun
    constants.mu_central = 132712439940;  % std. grav. parameter Sun  [km3 s-2]
    
    % standard body numbers
    % proper indices for arrays: Sun, Mercury, Venus, Earth, Mars, Jupiter,
    % Saturn, Uranus, Neptune, Moon, Ceres, Pallas, Vesta, Pluto
    model.numbers = 1:14;
    model.standard_bodies = 14;
    
    % model specific names
    model.names = {...
        'Sun';'Mercury';'Venus';'Earth';'Mars';'Jupiter';'Saturn';'Uranus';'Neptune';
        'Moon';'Ceres';'Pallas';'Vesta';'Pluto'};
    
    % mean radii ( http://ssd.jpl.nasa.gov, (volumetric) mean radius, [km] )
    Radii = [...Radii is re-used later on
        1.392e06; 2439.70; 6051.80; 6371.01; 3389.90; 69911.00; 58232.00; 25362.00; 24624.00;
        1737.53; 487.3; 580; 578; 1133];
    model.mean_Radii = Radii;
    
    % Minimal allowable altitudes
    model.min_alts = [...
        4*Radii(1); 200; 284; 306; 257; 42895; 20612; 81533; 4482;
        100; 0; 0; 0; 0];
    
    % default minimum transfer times
    % TODO: this is a rather difficult thing to default...Eventually these
    % numbers should be based on Hohmann transfer times or something, taking 
    % the next planet in the sequence...now we just use the same settings 
    % for all inner and outer planets...
    inner_LB = 25;   inner_UB = 500;  % days
    outer_LB = 350;  outer_UB = 2500; % days
    CBF_LB   = 25;   CBF_UB   = 750;  % days
    % default minimum transfer times
    model.TOF_LB = [CBF_LB;inner_LB;inner_LB;inner_LB;inner_LB;outer_LB;outer_LB;outer_LB;outer_LB;
        inner_LB;inner_LB;inner_LB;inner_LB;outer_LB];        
    % default maximum transfer times
    model.TOF_UB = [CBF_UB;inner_UB;inner_UB;inner_UB;inner_UB;outer_UB;outer_UB;outer_UB;outer_UB;
        inner_UB;inner_UB;inner_UB;inner_UB;outer_UB]; 
    
    % Is the body unmovable (convenient for mission targets that have a
    % fixed location w.r.t. the central body)
    model.static{1} = ...
        [ true; false; false; false; false; false; false; false; false;
        false; false; false; false; false]; 
    
    % if so, what are the coordinates?
    model.static{2} = {... % Sun is always at [0,0,0, 0,0,0]
        [0,0,0,0,0,0]; [];[];[];[];[];[];[];[];[];[];[];[];[]};
    
    % can the body be used as a launch site?
    model.departureable = ...
        [false; true; true;  true;  true;  true;  true;  true;  true;
        false; false; false; false; false];
    
    % Can the body be used for GAMs?
    model.GAMable = ...
        [ true;  true;  true;  true;  true;  true;  true;  true;  true;
        false; false; false; false; false];
    
    % Can the body be used for Aerogravity assists?
    model.aerogravable = ...
        [ false;  false;  true;  true;  true;  true;  true;  true;  true;
        false; false; false; false; false];
    
    % Can the body be a mission target?
    model.targetable = ...
        [ false;  true;  true;  true;  true;  true;  true;  true;  true;
        true; true; true; true; true];
    
end