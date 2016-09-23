function model = user_Solar_system_parameters(model)
% LOAD MISSION-SPECIFIC DATA
    
% location in ecliptic lon/lat
lambda = 75.4 *pi/180;
beta   = -7.5 *pi/180;

% XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx
% THESIS Stuff
% XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx

% stagnation point
one       = [lambda, beta];

% 6° perturbation
six_degrees = 6 *pi/180;
two       = [lambda, beta + six_degrees];
three     = [lambda - six_degrees/sqrt(2), beta + six_degrees/sqrt(2)];
four      = [lambda - six_degrees, beta];
five      = [lambda - six_degrees/sqrt(2), beta - six_degrees/sqrt(2)];
six       = [lambda, beta - six_degrees];
seven     = [lambda + six_degrees/sqrt(2), beta - six_degrees/sqrt(2)];
eight     = [lambda + six_degrees, beta];
nine      = [lambda + six_degrees/sqrt(2), beta + six_degrees/sqrt(2)];

% 12° perturbation
twelve_degrees = 12 *pi/180;
ten       = [lambda, beta + twelve_degrees];
eleven    = [lambda - twelve_degrees/sqrt(2), beta + twelve_degrees/sqrt(2)];
twelve    = [lambda - twelve_degrees, beta];
thirteen  = [lambda - twelve_degrees/sqrt(2), beta - twelve_degrees/sqrt(2)];
fourteen  = [lambda, beta - twelve_degrees];
fifteen   = [lambda + twelve_degrees/sqrt(2), beta - twelve_degrees/sqrt(2)];
sixteen   = [lambda + twelve_degrees, beta];
seventeen = [lambda + twelve_degrees/sqrt(2), beta + twelve_degrees/sqrt(2)];

% current test point
current = eleven;
% NOTE: Z = 0!!

% current test location
RSBS = 200*1.495978706910000e+008; 
X = RSBS*cos(current(2))*cos(current(1));
Y = RSBS*cos(current(2))*sin(current(1));
Z = 0;%RSBS*sin(current(2));
RSBS = [X,Y,Z];

% XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx
% THESIS Stuff
% XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx

    % modify appropriate parameters in [model]
    model.names         = [model.names; 'Solar bow shock']; % Solar bow shock (name)
    model.numbers       = (1:(numel(model.names)+1)).';     % Solar bow shock (reference)    
    model.static{1}     = [model.static{1}; true];          % SBS is always in the same place    
    model.static{2}     = [model.static{2}; 
                          [RSBS,...                         % (Fixed) statevector of SBS
                           0,0,0]];                          
    model.TOF_LB        = [model.TOF_LB; 1500];             % SBS upper and lower times of flight
    model.TOF_UB        = [model.TOF_UB; 5000];
    model.departureable = [model.departureable; false];     % SBS can NOT be used as a launch site
    model.GAMable       = [model.GAMable; false];           % SBS can NOT be used for swingbys
    model.aerogravable  = [model.aerogravable; false];      % SBS can NOT be used for aerograv-assists    
    model.targetable    = [model.targetable; true];         % SBS CAN be used as a mission target       

end 
