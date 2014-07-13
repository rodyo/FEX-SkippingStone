function varargout = modify_settings(funfcn, varargin)
    
    %% Initialize
    
    % get globals
    global MainWin launch_tab sequence_tab arrival_tab optimization_tab 
    
    % get all appdata
    if ishandle(MainWin)
        settings    = getappdata(MainWin, 'settings'   );
        handles     = getappdata(MainWin, 'handles'    );
        environment = getappdata(MainWin, 'environment');
        model       = getappdata(MainWin, 'model'      );
        constants   = getappdata(MainWin, 'constants'  );
        calculation = getappdata(MainWin, 'calculation');
    end
    
    % run appropriate Callback function         
    [ig,funfcn] = evalc(['@',funfcn]);%#ok               % (convert to function handle)
    [varargout{1:nargout}] = feval(funfcn, varargin{:}); % STR2FUNC() doesn't work in this context    
    
    % write all appdata back to MainWin
    if ishandle(MainWin)
        setappdata(MainWin, 'settings'   , settings   );
        setappdata(MainWin, 'handles'    , handles    );
        setappdata(MainWin, 'environment', environment);
        setappdata(MainWin, 'model'      , model      );
        setappdata(MainWin, 'constants'  , constants  );
        setappdata(MainWin, 'calculation', calculation);
    end
    
    %% Default values
    
    function [environment, model, constants, calculation, settings] = ...
            default_values(rootdir)%#ok
    % Load default values. This function is called on program
    % startup, or when the user selects the "reset to default"
    % menu option.
        
        %% environment
        
        % THE PROGRAM'S NAME
        environment.program_name = 'SkippingStone';
        
        % AUTHORS & CONTRIBUTORS
        % more authors can easily be added this way
        % just copy-paste this below, and fill in your own info:
        environment.authors = {};
        environment.authors = [environment.authors;{...
            'Rody P.S. Oldenhuis',...                       % name
            'oldenhuis@gmail.com',...                       % contact info
            'Delft University of Technology, Holland,',...  % affiliation field 1
            'Partly performed at ISAS/JAXA, Japan.',...     % affiliation field 2
            ' '}];                                          % affiliation field 3        
        % same thing for "contributors":
        environment.contributors = {...
            'John d`Errico',...
            'Yi Cao',...
            'Stephen Morris',...  
            'Erik Johnson',...
            'ioxv.4623 (all MATLAB FEX)',...
            'Dr. Dario Izzo (ESA/ACT)',...            
            };

        % determine whether Octave or MATLAB is running
        % (trick from author 'ioxv.4623' on matlab file exchange)
        environment.UI  = 'Unknown';
        LIC = license('inuse');
        for elem = 1:numel(LIC)
            envStr = LIC(elem).feature;
            if strcmpi(envStr,'Octave')
                environment.UI = 'Octave';
                break
            end
            if strcmpi(envStr,'Matlab')
                environment.UI = 'Matlab';
                break
            end
        end

        % pathing: set rootdir, datadir, previous path, and set proper path
        ep.rootdir     = rootdir;
        ep.datadir     = [rootdir, filesep, 'data'];
        ep.prevpath    = path;
        ep.MP_filename = [ep.datadir, filesep, ...
            'asteroids', filesep, 'MPCORB.DAT'];
        addpath(genpath(rootdir));
        environment.pathing = ep;

        % background color for windows
        environment.colors.window_bgcolor = get(0, 'defaultuicontrolbackgroundcolor');

        % NOTE - default BG-color for edit boxes is WHITE on Windows-systems.
        % and equal to the background color on Unix systems. Select the
        % appropriate value by using ISPC (see the help of ISPC or COMPUTER)
%         if ispc
            environment.colors.edit_bgcolor = 'White';
%         else
%             environment.colors.edit_bgcolor = ...
%                 get(0,'defaultUicontrolBackgroundColor');
%         end
        
        % see if the optimization toolbox is available
        if exist('fmincon', 'file') == 2
            environment.optim_toolbox_available = true;
        else
            environment.optim_toolbox_available = false;
        end
        
        %% plugins
        
        % load costfunctions (remove subdirs and extentions)
        costfuns = dir([rootdir filesep 'plugins' filesep 'costfunctions']);
        names = {costfuns.name};   dirs = [costfuns.isdir];
        costfuns = names(~dirs);   costfuns = cellfun(@(x)x(1:end-2),costfuns,'uniformoutput',false);
        
        % get info on each costfunction
        for i = 1:numel(costfuns)
            % do it carefully. If something goes wrong, make a report
            try
                % evaluate
                [dummy,costfun_info] = evalc(costfuns{i});%#ok
                % check the info
                if any(~isfield(costfun_info, ...
                        {'name', 'function_handle', 'description', 'axis_label'}))
                    error(' '); %<--- force the try-catch to fail
                end
                % insert in output argument.
                environment.plugin_info.costfuns(i) = costfun_info;
                
            catch %#ok
                errordlg(['Error detected in costfunction plugin ''', costfuns{i},...
                    ''': plugins must return ', 'a structure with the fields ''name'' ',...
                    '''function_handle'', ''description'' and ''axis_label'' when ',...
                    'called without input arguments.'], 'Plugin implementation error');
            end
        end
        
        % load post-processors (remove subdirs and extentions)
        postprocessors = dir([rootdir filesep 'plugins'  filesep 'postprocessors']);
        names = {postprocessors.name};  dirs = [postprocessors.isdir];
        postprocessors = names(~dirs);
        postprocessors = cellfun(@(x)x(1:end-2),postprocessors,'uniformoutput',false);
        
        % get info on each postprocessor
        for i = 1:numel(postprocessors)            
            % do it carefully. If something goes wrong, make a report
            try
                % evaluate
                [dummy, pp_info] = evalc(postprocessors{i});%#ok
                % check the info & GUI commands
                if any(~isfield(pp_info, {'name', 'function_handle', ...
                        'GUI_function_handle' 'plot_function_handle'}))
                    error(' '); %<--- force the try-catch to fail
                end
                % insert in output argument
                environment.plugin_info.postprocessors(i) = pp_info;
                                
            catch %#ok
                errordlg(['Error detected in post-processor plugin ''', postprocessors{i},...
                    ''': plugins must return a structure with the fields ''name'', ',...
                    '''function_handle'', ''GUI_function_handle'' and ',...
                    '''plot_function_handle'' when called without input arguments.'], ...
                    'Plugin implementation error');
            end
        end
        
        %% model & constants

        % load default model & constants
        [model, constants] = Solar_system_parameters([],[]);
        
        % load mission-specific data for the Solar system
        model = user_Solar_system_parameters(model);

        % Astrodynamical constants from (From NASA/JPL Horizons)
        % and some convenient constants
        constants.G  = 6.67300e-20;   % gravitational constant    [km+3 kg-1 s-2]
        constants.AU = 149597870.691; % Astronomical Unit         [km]
        constants.g0 = 9.80665/1e3;   % grav. acc. at sealevel    [km s-2]
        % (note the conversion to km)
        % seconds, minutes, hours, days, weeks, months, years
        constants.timeunits = {1; 60; 3600; 86400; 7*86400; 30.471*86400; 365.25*86400};
        % keep this a Cell-array; it allows one to use DEAL()

        %% calculation
        
        % initialize first-order results structure
        c1.Paretos.solutions       = [];        c1.Paretos.function_values = [];    
        c1.best.solution           = [];        c1.best.function_value     = [];
        c1.elapsed_time            = [];        c1.lambert_problems        = [];
        c1.ephemerides             = [];        c1.algorithm               = [];
        
        % initialize BATCH
        calculation.results.BATCH = [];
        
        % initialize second-order results structure               
        c2.elapsed_time  = [];        c2.best.solution = [];
        c2.integrations  = [];        c2.algorithm     = [];
        c2.ephemerides   = [];    
        
        % insert them into calculation
        calculation.results.first_order  = c1; clear c1;
        calculation.results.second_order = c2; clear c2;
        
        
        %% settings, Launch & Satellite Data
                
        % launch
        sl.launch_mass  = 3000;       % wetmass                        [kg]
        sl.launch_mass_margin = 0;    % wetmass margin                 [%]
        sl.adapter_mass = 0;          % adapter mass                   [kg]
        sl.max_C3 = 225;              % max. C3 at launch              [km2/s2]
            % when launchers get implemented, this should become a vector,
            % the N-th value of which corresponds to the Launcher        
        sl.payload_mass = 50;         % minimum drymass                [kg]
        sl.payload_mass_margin = 0;   % drymass margin                 [%]
        % insert these settings
        settings.launch = sl; clear sl;
        % power
        sp.selected = [1 0];           % which one's selected
        % power (SEP)
        sp.SEP.power_1AU = 10000;      % power at 1 AU                  [W]
        sp.SEP.deterioration_rate = 0.9883; % deterioration rate        [-]
            % deterioration is computed as P(t) = P0 * rate^t, [t] in years
        sp.SEP.jettisonable = false;   % panels jettisonable?           [-]
        sp.SEP.panel_mass = 50;        % mass of panels                 [kg]
        % power (NEP)
        sp.NEP.amount = 1;             % amount of Nuclear batteries
        sp.NEP.power_per_battery = 150;% power per battery              [W]
            % when different battery types get implemented, this should 
            % become a vector, the N-th value of which corresponds to 
            % the power produced by that model
        sp.NEP.deterioration_rate = -0.9; % battery deterioration rate  [-]
            %??? not yet implemented:  P(t) = P0 * exp(rate * t)
        % power (general)
        sp.instrument_power = 0;       % minimal power for instruments  [W]
        % insert these settings
        settings.power = sp; clear sp;
        % thrust & engine data
        sp.selected = [1 0 0];    % Which one's selected? 
        % thrust & engine data, high thrust        
        sp.high_thrust.Isp = 300; % Isp for chemical thruster      [s] 
        % thrust & engine data, ion engine
        sp.ion_engine.Isp = 3000;   % Isp for ion engine           [s] 
        sp.ion_engine.maxT = 0.135; % maximum thrustlevel          [N]
        sp.ion_engine.amount = 1;   % number of engines            [-]
        sp.ion_engine.min_power = 0.649e3; % min. power level      [W]
        sp.ion_engine.duty_cycle = 0.9; % duty cycle               [%]
            % ??? not yet implemented
        settings.propulsion.ion_engine.efficiency = 0.8; % efficiency of the engine [%]
            % ??? not yet implemented
        % Solar sail
        sp.Solar_sail.surface_area = 250;  % sail surface area      [m2]
        sp.Solar_sail.reflectivity = 0.98; % sail reflectivity      [%]
            % ??? not yet implemented
        % different Isp for capture
        sp.different_Isp.check = false; % use different Isp?
        sp.different_Isp.Isp = 300;     % default different Isp     [s]        
        % insert these settings
        settings.propulsion = sp; clear sp;
        

        %% settings, Sequence
        
        % initialize
        environment.numGAMpanels = 4;            % number of swingby-panels
        one = ones(environment.numGAMpanels, 1); % facilitates assignments        
        Earth = 3;                               % index to Earth;        
        % departure
        sd.body = Earth;                  % defaults to Earth
        sd.launch_window.day(1)   = 1;    % launchwindow, start day
        sd.launch_window.month(1) = 1;    % launchwindow, start month
        sd.launch_window.year(1)  = 66;   % launchwindow, start year (N-th entry)
        sd.launch_window.day(2)   = 1;    % launchwindow, end day
        sd.launch_window.month(2) = 1;    % launchwindow, end month
        sd.launch_window.year(2)  = 76;   % launchwindow, end year (N-th entry)
        % insert these settings
        settings.departure = sd; clear sd;
        % model
        settings.model = [true; false; false; false; false]; % selected model
            % format: [Solar Jovian Julian MPCORB user]
        % swingbys
        % default values for upper and lower bounds on TOF
        sG.body = 1*one;           % no bodies selected    
        sG.TOF_LB = [0; model.TOF_LB]; % time of flight, lower bounds (include zero for 'none' option)
        sG.TOF_UB = [0; model.TOF_UB]; % time of flight, upper bounds (include zero for 'none' option)
        sG.TOF_LB = sG.TOF_LB(:, one); % replicate for every GAM-panel
        sG.TOF_UB = sG.TOF_UB(:, one); % replicate for every GAM-panel
        sG.min_altitude = repmat(model.min_alts, 1, environment.numGAMpanels); 
                                               % minimum altitudes (use model default)         
        sG.short_or_longway = 3*one; % optimize short or long way         
        sG.max_DV = 5*one;         % max. Delta-V (100 is unlimited)
        sG.complete_revolutions  = 1*one; % use zero complete revolutions
        sG.LoverD = 2*one;           % L / D parameter [-]
        sG.type   = repmat({'powered'}, numel(one),1); % use powered swingby's 
        sG.jettison = false*one;     % don't jettison the panels anywhere        
        % swingby constraints
        sG.constraints.max_tof = 9131;           % maximum mission duration [days]
        sG.constraints.max_DV  = 10;              % maximum total DeltaV     [km/s]
        sG.constraints.min_Solar_distance = 0.3;  % minimum Solar distance   [AU]
        sG.constraints.C3LoverD_tolerance = 5;    % C3 / L/D tolerance       [%]
        % insert these settings
        settings.GAM = sG; clear sG;
        % target
        st.body = 1;     % first body selected
        st.TOF_LB = model.TOF_LB(model.targetable); % time of flight, lower bounds
        st.TOF_UB = model.TOF_UB(model.targetable); % time of flight, upper bounds 
        st.complete_revolutions = 1; % use zero complete revolutions
        st.short_or_longway = 3;     % optimize
        % insert these settings
        settings.target = st; clear st;
        % batch optimization 
        settings.BATCH.check = false; % default is OFF
        
        
        %% settings, Arrival data & postprocessing
        
        % arrival
        sa.type = 1;            % arrival type is flyby         
        sa.user_cost_function = false; % use user-defined cost function at arrival
        
        % arrival constraints
        sa.constraints.max_C3 = inf; % max C3 at arrival        
        % ??? TODO: these should depend on target
        sa.constraints.apocenter_altitude  = 250;
        sa.constraints.pericenter_altitude = 250;
        sa.constraints.inclination         = 0;
        
        sa.gravity_loss = 0;    % gravity loss at arrival      [%]
        sa.mass_released = 0;   % mass released before arrival [kg]
        
        % post-processing panel
        spp.check = false;
        spp.post_processor = 1; % (none) by default        
        
        % "Find MP's close to the trajectory" post-processor
        spp.optimize_distance.check        = true;
        spp.threshold.check(1)             = false;
        spp.constant_threshold.distance    = 0.01;
        spp.threshold.check(2)             = true;
        spp.variable_threshold.distance(1) = 0.01;
        spp.variable_threshold.distance(2) = 0.09;
        spp.variable_threshold.distance(3) = 30;
        
        % insert these settings 
        settings.postprocessing = spp; clear spp;
        settings.arrival        = sa;  clear sa;
      
        %% settings, Optimization 
        
        so.global.optimizer = [1,0,0];       % GODLIKE
        so.global.low_thrust_approximation = [1,0]; % Exponential Sinusoids        
        % standards for OPTIMIZE()
        oo = optimset('maxfunevals', 1e6, 'maxiter', 1e3); 
        oo.TolCon = 1e-6;  oo.popsize = [];  % non-standard options   
        so.global.optimizer_settings{3} = oo;% standard OPTIMIZE (fminsearch) options
        so.global.optimizer_settings{2} = oo;% standard OPTIMIZE (fminblfgs) options
        so.global.optimizer_settings{1} = set_options; % standard GODLIKE options
        so.ephemerides = [1,0,0];            % JPL / DE405 ephemerides
        so.objectives.max_mass = true;       % optimize for maximum payload mass
        so.objectives.min_tof = false;       % don't optimize for min. time-of-flight
        so.objectives.other.use = false;     % don't optimize for custom objectives
        so.objectives.other.function_handle = [];   % placeholder for the function to use
        so.objectives.other.name     = '';   % placeholder for the function name        
        so.objectives.other.axis_label = ''; % placeholder for the axis label
        if environment.optim_toolbox_available
            so.local.optimizer = [1,0,0];    % FMINCON is possible
        else
            so.local.optimizer = [0,0,1];    % otherwise, Quasi-Newton L-BFGS
        end
        so.local.integrator = [0,1,0];       % RKN8(6) integrator
        so.local.method     = [1,0];         % patched micro-conics        
        % insert these settings into settings
        settings.optimize = so;
        
        %% settings, Output
        
        % no settings for output
        
    end % default values

    %% Change all settings
    
    % set all values according to current settings
    function change_all_settings%#ok
        
        %% Launch & Satellite Data
        
        % some abbreviations
        lt  = handles.tab(launch_tab);
        sl  = settings.launch;
        sp  = settings.power;
        spp = settings.propulsion;
        
        % launch
        set(lt.LaunchMass          , 'string', sl.launch_mass       );
        set(lt.LaunchMassMargin    , 'string', sl.launch_mass_margin);
        set(lt.AdapterMass         , 'string', sl.adapter_mass      );
        set(lt.MaxC3Launch         , 'string', sl.max_C3            );  
            %??? this should be changed when implementing launchers      
        set(lt.MinPayloadMass      , 'string', sl.payload_mass      );
        set(lt.MinPayloadMassMargin, 'string', sl.payload_mass_margin);               
        % power
        for ii = 1:2, set(lt.power_radio(ii), 'value', sp.selected(ii)); end
        % power (SEP)
        set(lt.SEP(2)     , 'string', sp.SEP.power_1AU         );
        set(lt.SEP(5)     , 'string', sp.SEP.deterioration_rate);
        set(lt.SEP(7)     , 'value',  sp.SEP.jettisonable      );
        set(lt.Jettison(1), 'string', sp.SEP.panel_mass        );
        % power (NEP)
        set(lt.NEP(2),  'string', sp.NEP.amount           );
        set(lt.NEP(5),  'string', sp.NEP.power_per_battery);
            %??? this should be changed when implementing different battery models      
        %set(???, sp.NEP.deterioration_rate);            
        % power (general)
        set(lt.MinInstrumentPower, 'string', sp.instrument_power);
        % thrust & engine data        
        for ii = 1:3, set(lt.propulsion(ii), 'value', spp.selected(ii)); end
        if get(lt.propulsion(1), 'value')
            set(lt.High(2), 'string', spp.high_thrust.Isp);
        elseif get(propulsion(2), 'value')
            set(lt.High(2), 'string', spp.ion_engine.Isp );
        end
        set(lt.Ion(2) , 'string', spp.ion_engine.maxT        );
        set(lt.Ion(5) , 'string', spp.ion_engine.amount      );
        set(lt.Ion(8) , 'string', spp.ion_engine.min_power   );         
        %set(???, 'string', spp.ion_engine.duty_cycle );
        %set(???, 'string', spp.ion_engine.efficiency );
        set(lt.Sail(2), 'string', spp.Solar_sail.surface_area);
        %set(???, 'string', spp.Solar_sail.reflectivity);
        set(lt.Isp(1) , 'value' , spp.different_Isp.check    );
        set(lt.Isp(3) , 'string', spp.different_Isp.Isp      ); 
        
        % just to be sure
        clear lt sl sp spp;
        
        %% Sequence
        
        % some abbreviations
        st = handles.tab(sequence_tab);    sG = settings.GAM;
        sT = settings.target;              sd = settings.departure;

        % departure        
        for ii = 1:2
            set(st.launch_window.year(ii) , 'value', sd.launch_window.year(ii));
            set(st.launch_window.month(ii), 'value', sd.launch_window.month(ii)      );
            set(st.launch_window.day(ii)  , 'value', sd.launch_window.day(ii)        );
        end
        set(st.DepartureBody, 'value', sd.body);         
        % model
        for ii = 1:5, set(st.selected_model(ii), 'value',  settings.model(ii)); end
        % swingbys
        for ii = 1:environment.numGAMpanels
            set(st.GAM.body(ii)     , 'value' , sG.body(ii)            );
            set(st.GAM.TOF_LB(ii)   , 'string', sG.TOF_LB(sG.body(ii),ii)); 
            set(st.GAM.TOF_UB(ii)   , 'string', sG.TOF_UB(sG.body(ii),ii)); 
            set(st.GAM.shortlong(ii), 'value' , sG.short_or_longway(ii));
            set(st.GAM.maxDeltaV(ii), 'string', sG.max_DV(ii)          );
            set(st.GAM.CompleteRevolutions(ii), 'value', sG.complete_revolutions(ii));
            set(st.GAM.LoverD(ii)   , 'string', sG.LoverD(ii)          ); 
            current_GAMtypes = get(st.GAM.type(ii), 'string');
            set(st.GAM.type(ii)     , 'value' , find(strcmpi(sG.type{ii}, current_GAMtypes)));            
            set(st.GAM.jettison(ii) , 'value' , sG.jettison(ii)        );
            set(st.GAM.minalt(ii)   , 'string', sG.min_altitude(ii)    );
        end  
        % swingby constraints
        set(st.MaxTotalTOF      , 'string', sG.constraints.max_tof           );
        set(st.MaxTotalDeltaV   , 'string', sG.constraints.max_DV            );
        set(st.MinDistToSun     , 'string', sG.constraints.min_Solar_distance);
        set(st.C3LoverDTolerance, 'string', sG.constraints.C3LoverD_tolerance);
        % target
        set(st.target.body  , 'value' , sT.body  );
        set(st.target.TOF_LB, 'string', sT.TOF_LB(sT.body));
        set(st.target.TOF_UB, 'string', sT.TOF_UB(sT.body));
        set(st.target.CompleteRevolutions, 'value', sT.complete_revolutions);
        set(st.target.shortlong, 'value', sT.short_or_longway);

        % batch optimization 
        set(st.BatchOptimizeCheck, 'value', settings.BATCH.check);
        
        % just to be sure
        clear sT st sG sd;
        
        %% Arrival data & postprocessing

        % some abbreviations
        at  = handles.tab(arrival_tab);        sa = settings.arrival;
        spp = settings.postprocessing;
        
        % arrival
        set(at.ArrivalType           , 'value' , sa.type              );        
        set(at.usr_costfun_arrival(1), 'value' , sa.user_cost_function);
        
        % arrival constraints
        set(at.MaxArrivalC3        , 'string', sa.constraints.max_C3);
        
        set(at.ApocenterAltitude   , 'string', sa.constraints.apocenter_altitude );
        set(at.PericenterAltitude  , 'string', sa.constraints.pericenter_altitude);
        set(at.Inclination         , 'string', sa.constraints.inclination);
        
        set(at.ArrivalGravityLoss  , 'string', sa.gravity_loss      );
        set(at.InsertionMassRelease, 'string', sa.mass_released     );
                
        % post-processing
        set(at.post_processing.check, 'value', spp.check);
        set(at.post_processing.post_processor(2), 'value', spp.post_processor);
        
        % just to be sure
        clear at spp sa;
        
        %% Optimization
        
        % some abbreviations
        ot = handles.tab(optimization_tab);
        so = settings.optimize;
        
        % data
        for ii = 1:3, set(ot.global_optimizer(ii), 'value', so.global.optimizer(ii)); end
        for ii = 1:2, set(ot.LowThrustApproximation(ii), 'value', so.global.low_thrust_approximation(ii)); end
        for ii = 1:3, set(ot.Ephemerides(ii), 'value', so.ephemerides(ii)); end
        set(ot.MaxPayloadObjective  , 'value', so.objectives.max_mass);
        set(ot.MinTOFObjective      , 'value', so.objectives.min_tof );
        set(ot.OtherObjectives.check, 'value', so.objectives.other.use);
        % re-check if optimization toolbox is available
        environment.optim_toolbox_available = (exist('fmincon', 'file') == 2);                   
        if environment.optim_toolbox_available && so.local.optimizer(1)
            set(ot.local_optimizer(1), 'value', 1);
        elseif settings.optimize.local.optimizer
            set(ot.local_optimizer(2), 'value', 1);            
        else
            for ii = 1:3, set(ot.local_optimizer(ii), 'value', so.local.optimizer(ii)); end
        end
        for ii = 1:3, set(ot.integrator(ii), 'value', so.local.integrator(ii)); end
        for ii = 1:2, set(ot.low_thrust_optimization_method(ii), 'value', so.local.method(ii)); end
        
        % just to be sure
        clear ot so;

        %% Output
        
        % nothing to load
  
        % environment might have changed, so insert it     
        setappdata(MainWin, 'environment', environment);          
        
    end % change all settings
    
    %% Change single setting
    
    % change a single setting (general callback for all controls)
    function stdargout = change_single_setting(setting, handle, varargin)%#ok
        
        %% initialize
        
        % this function has to assign at least one value
        % (allows multiple callbacks per UICONTROL())
        stdargout = [];   
        % rename calling object        
        object = varargin{1};        
        % get its type and its new setting
        try
            % editboxes, checkboxes and popupmenus  
            type = get(object, 'style');            
        catch %#ok
            % buttongroups
            type = get(object, 'type');
        end        

        %% edit boxes
        
        if strcmpi(type, 'edit')
            % extract new setting
            new_setting = get(object, 'string');           
            % check inserted value
            good = check_value(str2double(new_setting));
            % insert or reject the new setting
            if good
                evalc(['settings.', setting, '=', new_setting]);                
            else
                evalc(['set(object, ''string'', settings.', setting,')']);
                return
            end   

        %% check boxes, popups
        
        elseif any(strcmpi(type, {'checkbox';'popupmenu'}))
            % extract new setting
            new_setting = get(object, 'value');
            % insert the new setting
            evalc(['settings.', setting, '=', num2str(new_setting)]);

        %% buttongroups
        
        elseif strcmpi(type, 'uipanel')
            % find which one's selected
            selected_radio = get(object, 'selectedobject');%#ok
            [dontcare, which_one] = evalc([handle, ' == selected_radio']);%#ok
            % insert new data
            evalc(['settings.', setting, '=[', num2str(which_one),']']);
        end     
    end  % change single setting
    
end % all settings-related functions

%% Check input from editboxes

function good = check_value(value)
    % optimistic default
    good = true;            
    % only numbers are allowed
    if ~isfinite(value) % STR2DOUBLE() returns NaN for non-numeric input
        uiwait(errordlg({'Input can only contain numbers 0-9, periods (.),';
            'exponents (e) and signs (+/-).'}, 'Only numbers allowed')); uiresume;
        good = false; return;
    end        
    % only REAL numbers are allowed
    if ~isreal(value) % I've had a "+i" typo once...don't ask
        uiwait(errordlg('Value must be real.', 'Complex value received')); uiresume;
        good = false; return;
    end        
    % only POSITIVE numbers are allowed
    if (value < 0)
        uiwait(errordlg('Value must be positive.', 'Negative value received')); uiresume;
        good = false; return;
    end        
end
