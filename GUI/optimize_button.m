%% Callback function for the optimize! button

function optimize_button(varargin) %#ok<VANUS>

    %% Initialize

    % get globals
    global MainWin sequence_tab output_tab algorithms_tab

    % extract relevant data
    settings    = getappdata(MainWin, 'settings'   );
    calculation = getappdata(MainWin, 'calculation');
    handles     = getappdata(MainWin, 'handles'    );
    environment = getappdata(MainWin, 'environment');

    %% Perform basic checks

    % wetmass lower than drymass
    if (settings.launch.launch_mass < settings.launch.payload_mass)
        uiwait(errordlg('Launch mass can not be lower than minimum drymass.',...
            'Incorrect parameter values', 'modal')); uiresume % <- prevents interaction with
                                                              %    GUI until the error is OK'ed
        return;
    end

    % no swingby-bodies selected
    swingby_bodies = get(handles.tab(sequence_tab).GAM.body, 'value');
    if sum([swingby_bodies{:}] - 1) < 1
        errordlg({'You have not selected any swinbgby bodies.';
                 [environment.program_name, ' can only optimize MGA-problems.']},...
                 'No swinbgby bodies selected')
        return
    end

    % convert launch window dates to MJD
    initial_year = get(handles.tab(sequence_tab).launch_window.year(1), 'string');
    initial_year = str2double(initial_year(1, 1:4));
    launch_window_days_past_J2000 = date2days(settings.departure.launch_window.year + initial_year - 1,...
                                              settings.departure.launch_window.month,...
                                              settings.departure.launch_window.day,...
                                              0, 0, 0);

    % launch window end-date preceeds start-date
    if diff(launch_window_days_past_J2000) < 0
        uiwait(warndlg({'Launch window end date preceeds start date;';
                       'values will simple be swapped.'},...
                       'Launch window dates reversed',...
                       'modal')), uiresume; % <- blocks execution until OK is pressed
        % swap values
        temp = launch_window_days_past_J2000(1);
        launch_window_days_past_J2000(1) = launch_window_days_past_J2000(2);
        launch_window_days_past_J2000(2) = temp;
    end

    % insert MJD into settings
    settings.launch_window_days_past_J2000 = launch_window_days_past_J2000;

    % store it in AppData
    setappdata(MainWin, 'settings', settings);

    % First (global) or second (local) order optimization?
    if ~isempty(calculation.results.first_order.best.solution)
        button_pressed = questdlg(['A result from a previous global optimization still resides in memory. ',...
                                  'You may choose to use this result as the initial value for a high-accuracy local search, ',...
                                  'or perform a new medium-accuracy global search with the current settings.'],...
                                  'Choose optimiation type',...
                                  'Perform new global search', ...
                                  'Perform high-accuracy local search',...
                                  'Perform high-accuracy local search');

        switch lower(button_pressed)

            case 'perform new global search'
                % clearing this variable will cause MGA to re-run the
                % global search; this is the only criterion to enter
                % that part of the function
                calculation.results.first_order.best.solution = [];
                setappdata(MainWin, 'calculation', calculation);

                toggle_panel(handles.tab(algorithms_tab).second_order_panel, 'disable');

            case 'perform high-accuracy local search'
                % do nothing; a non-empty [...best.solution] will cause
                % MGA to use the best solution for the second-order
                % search

            otherwise
                return;
        end
    end

    %% Start optimization

    % disable all controls in the Main window
    % (cancel button is enabled AFTER the models have loaded)
    [objects, states] = callbacks('disable_all');

    % Make 100% sure that the Cancel-button state is neutral
    set(handles.CancelButton, ...
        'UserData', false,...
        'enable'  , 'on');

    % set MainWindow flag to dirty (unsaved optimization)
    set(MainWin, 'UserData', 'dirty');

    % run the optimization
    try
        % set the progress bar
        progress_bar(0, 'Initializing optimization...')

        % BATCH-optimizations
        if settings.BATCH.check
            %{
            %TODO: implement BATCH optimizations

            % For BATCH optimizations, we first have to select the
            % different parameters according to each sequence. Also, the
            % optimization results need to be saved differently.

            % loop through all sequences
            for ii = 1:settings.BATCH.num_combinations

                % current sequence
                seq = settings.BATCH.combinations{ii};
                % select (default) times of flight
                %??? TODO

                % use "best-of-X" optimizations
                result = cell(setting.BATCH.best_of,1);
                for jj = 1:setting.BATCH.best_of
                    [result{jj}, successful] = MGA('batch'); end %#ok<NASGU>

                % find the best result & insert in calculation
                %??? TODO
                calculation.results.BATCH(ii) = best_result;

            end % loop through all BATCH-sequences

            % select the VERY best result
            %??? TODO

            % there's bound to be ONE successful
            successful = true;

            %}

        % single optimziations
        else
            % Call the main function
            do_again = true;
            while do_again
                [calculation, successful, do_again] = MGA(); end
        end

    % optimization might fail for some un-debugged reason
    catch ME  %#ok<MUCTH>
        
        % no success
        successful = false;
        
        % form MATLAB-exception object
        mME = MException('MainWin:optimization_failed', ...
                         'Unhandled exception; this is a bug.');                     
        ME = addCause(ME, mME);
        
        % report error
        report = getReport(ME, 'extended', 'hyperlinks', 'off');
        
        % show warning
        uiwait(warndlg({'Something went wrong:'; ' ';  report}, ...
                        'Optimization failed!', ...
                        'modal')); 
        uiresume();
    end

    % reset all controls in the Main Window
    callbacks('reset_all', objects, states);

    %% Process results

    % check if optimization was cancelled
    if cancel_button_pressed()
        % first write results to appdata
        setappdata(MainWin, 'calculation', calculation);
        % adjust progress bar and show warning
        progress_bar(0, 'Operation aborted')
        uiwait(warndlg({'Optimization aborted!';
            'Best result(s) found thus far will be used.'},...
            'Optimization cancelled', ...
            'modal')); uiresume % <- blocks execution until OK is pressed
        % but DO plot results
        callbacks('showtab', output_tab); % switch to output tab
        generate_output('embedded');      % and display results
        % reset everything
        set(handles.CancelButton, 'UserData', 0);
        progress_bar('');

    % check if optimization failed
    elseif ~successful
        % reset progress bar
        progress_bar('')
        %??? TODO

    % otherwise, optimization completed succesfully
    else
        % first write results to appdata
        setappdata(MainWin, 'calculation', calculation);

        % enable local search panel on the Algorithms tab
        toggle_panel(handles.tab(algorithms_tab).second_order_panel, 'reset');

        progress_bar(1, 'Optimization completed sucessfully.')

        % finish off this optimization
        callbacks('showtab', output_tab); % switch to output window
        generate_output('embedded');      % and display results

        progress_bar('');
    end

end % optimize button

%% MGA (Multiple Gravity Assist) - wrapper function
%
% This function simply extracts all current settings from the Main
% Window and uses these to construct the correct input for all the
% external functions. Both the global and high-accuracy local
% optimization are performed here.
function varargout = MGA(batch_or_normal)

    %% Initialize

    % get handle to MainWin
    global MainWin

    % extract data
    settings    = getappdata(MainWin, 'settings'   );
    model       = getappdata(MainWin, 'model'      );
    environment = getappdata(MainWin, 'environment');
    calculation = getappdata(MainWin, 'calculation');
    constants   = getappdata(MainWin, 'constants'  );

    % Determine GAM-sequence (remove 'none')
    seq = settings.GAM.body - 1;
    first_none = find(seq == 0, 1);
    if ~isempty(first_none)
        seq = seq(1:first_none-1); end

    % amount of swingby's to perform
    num_swingbys = numel(seq);
    % append departure and target bodies
    % (+1 because the central body can NOT be used as departure or
    % target body, so it is NOT included in the corresponding lists)
    seq = [settings.departure.body+1; seq; settings.target.body+1];

    %% Load correct model and initialize ephemerides

    % Solar system model
    if settings.model(1)
        % initialize the ephemerides model and wrapper function
        model = initialize_ephemerides_generators(seq, model, environment, 'Solar_System');
        ephemerides_generation([],[], model, environment, settings);
        % load complete solar system model
        model = Solar_system_model(seq, ...
            settings.launch_window_days_past_J2000(1), model, environment, constants);
        % load planetary atmospheres model (when required)
        if any(strcmpi(settings.GAM.type(1:num_swingbys), 'aerograv'))
            model = planetary_atmospheres_model(seq, model);
        end

    % Jovian model
    elseif settings.model(2)
        % initialize the ephemerides model and wrapper function
        model = initialize_ephemerides_generators(seq, model, environment, 'Jovian_System');
        ephemerides_generation([],[], model, environment, settings);
        % load complete Jovian system model
        model = Jovian_system_model(model);
        %??? to be continued later

    % Julian model
    elseif settings.model(3)
        % initialize the ephemerides model and wrapper function
        model = initialize_ephemerides_generators(seq, model, environment, 'Julian_System');
        ephemerides_generation([],[], model, environment, settings);
        % load complete Jovian system model
        model = Julian_system_model(model);
        %??? to be continued later

    end

    % Minor planets
    if settings.model(4)
        % load the database
        [model, constants] = minor_planets_model(model, constants, environment);
        % FIRST save the model - user functions MIGHT contain errors (^_^), in
        % which case the ENTIRE model is lost; the program would need to be
        % restarted, and the MPCORB-database reloaded (annoying!!)
        setappdata(MainWin, 'model', model);
        % Now load the user-defined modifications
        model.MPs = user_MP_model(model.MPs);
        % and apply mission-specific pruning rules
        model.MPs = user_MP_pruning(model.MPs, constants);
    end

    % User database
    if settings.model(5)
        %??? TODO: to be implemented later
    end

    % Make 100.4% sure the full model is saved
    setappdata(MainWin, 'model', model);

    %% Parse current settings

    % set standard options for the patched conics procedure

    % WARNING: don't insert [model] directly -- without MP's it's no
    % problem, but when the MP's are included in the model (post-
    % processors/MP-costfunctions), more time is actually spent on
    % copying [model] to the ephemerides generator than on the actual
    % optimization. This is because the complete MP-model is over 1 GB
    % in size.
    %
    % NOTE: you are probably reading this because you ended up here
    % because of an error regarding some indexing operation. You will
    % get this error because you indlucded some swingby body that is
    % NOT part of the standard model you are using. If this is true,
    % you can correct this error by defining the following parameters
    % for your non-standard bodies (in the appropriate
    % /MISSION_SPECIFIC/USR_<model>_model.m file):
    %
    %   VescSOI         Vesc at the sphere of influence
    %   GMs             std. grav. parameter (=GM)
    %   mean_Radii      mean radius of your body
    %   min_altitude    minimum allowable altitude
    %
    % Of course, one (or more) of these may not be applicable for your
    % body. In that case, just define a value of [inf], [NaN] or 0. As
    % long as the position in the appropriate array is filled, the error
    % will not occur.
    pc_params = struct(...
        'scope'               , 'local',...
        'ephemerides_generator', @(body, time) ephemerides_generation(body, time),...
        'M0'                   , settings.launch.launch_mass,...
        'Me'                   , settings.launch.payload_mass,...
        'VescSOI'              , model.VescSOI(seq(2:end-1)),...
        'GMs'                  , model.GMs(seq(2:end-1)),...
        'muC'                  , constants.mu_central,...
        'mean_Radii'           , model.mean_Radii(seq(2:end-1)),...
        'max_total_DV'         , settings.GAM.constraints.max_DV,...
        'max_C3'               , settings.launch.max_C3,...
        'max_arrival_C3'       , settings.arrival.constraints.max_C3,...
        'C3LoverD_tolerance'   , settings.GAM.constraints.C3LoverD_tolerance,...
        'min_alts'             , settings.GAM.min_altitude(seq(2:end-1) + ...
                                 model.standard_bodies*(0:num_swingbys-1).'),...
        'types'                , {settings.GAM.type(1:num_swingbys)},...
        'mindist_to_central'   , settings.GAM.constraints.min_Solar_distance*constants.AU,...
        'initial_mass'         , settings.launch.launch_mass,...
        'max_GAM_DeltaV'       , settings.GAM.max_DV(1:num_swingbys),...
        'max_TOF'              , settings.GAM.constraints.max_tof);

    % modify PC_PARAMS according to intended type of optimization
    if isempty(calculation.results.first_order.best.solution)
        %% First-order (Global optimization)

        % determine if high or low thrust has been selected
        high_thrust = settings.propulsion.selected(1);
        % optimize Globally
        first_order = true;   second_order = false;

        % single- or multi-objective?
        single_objective = nnz(...
            [settings.optimize.objectives.max_mass
            settings.optimize.objectives.min_tof
            settings.optimize.objectives.other.use]) == 1;
        multi_objective = ~single_objective;

        % high-thrust
        if high_thrust
            pc_params.solution_type = 'ballistic';
            pc_params.Isp = settings.propulsion.high_thrust.Isp;

        % low-thrust
        else
            % ExpoSins
            if settings.optimize.global.low_thrust_approximation(1)
                pc_params.solution_type = 'ExpoSins';
                pc_params.Isp = settings.propulsion.ion_engine.Isp;
            % equinoctial elements
            elseif settings.optimize.global.low_thrust_approximation(2)
                pc_params.solution_type = 'Equinoctial elements';
                pc_params.Isp = settings.propulsion.ion_engine.Isp;
            end

            %??? Solar sail lambert targeter?

        end

    % Second-order (Local optimization)
    else
        %% Second-order (Local optimization)

        % determine if high or low thrust has been selected
        high_thrust = settings.propulsion.selected(1);
        % optimize Locally
        first_order = false;   second_order = true;

        % high-thrust
        if high_thrust
            pc_params.solution_type = 'ballistic-integration';
            pc_params.Isp = settings.propulsion.high_thrust.Isp;
            % first make sure the ENTIRE model is loaded
            models = {'Solar_System'; 'Jovian system'; 'Julian system'};
            model = initialize_ephemerides_generators(...
                1:model.standard_bodies, model, models{settings.model(1:3)});
            setappdata(MainWin, 'model', model);
            % selected type of integrator
            integrators = {...
                @(odefun,tspan,y0,options)      ode113(odefun, tspan, y0, options)
                @(odefun,tspan,y0,dy0,options)   rkn86(odefun, tspan, y0, dy0, options)
                @(odefun,tspan,y0,dy0,options) rkn1210(odefun, tspan, y0, dy0, options)
                @(odefun,tspan,y0,dy0,options) []   % TODO: ODEX2
                @(odefun,tspan,y0,dy0,options) []}; % TODO: GBS
            integrator = integrators{logical(settings.optimize.local.integrator)};
            % pass all data also to NBody()
            int_options = struct(...
                'GMs'        , model.GMs,...  % std. grav. parameters of all perturbers
                'integrator' , integrator,... % type of integrator
                'perturbers' , 1:model.standard_bodies,...       % ALL model-bodies are perturbers
                'ephemerides', pc_params.ephemerides_generator); % copy ephemerides generator
            % and set integrator
            pc_params.integrator = @(y0, tspan) NBody_Battin(y0, tspan, int_options);

        % low-thrust
        else
            % Sims & Flanagan

            % Collocation
        end

    end

    %% First-order optimization

    if first_order

        %% optimization & post-processing

        % initialize
        progress_bar(0, 'Optimizing...')                 % update progress bar
        pause(0.25), start_time = tic;                   % start timer
        PC_objective_function([],[], 'reset_stats');     % reset persistent variables
        objectives.which = [...                          % select which objectives are to be optimized
            settings.optimize.objectives.max_mass        % maximum mass objective
            settings.optimize.objectives.min_tof         % minimum time of flight objective
            settings.optimize.objectives.other.use];     % other objectives
        objectives.other = settings.optimize.objectives.other.function_handle;

        % objective function for global optimizer
        objective_function = @(X, params) PC_objective_function(...
            objectives, patched_conics(seq, X, params));

        % call single-objective optimizer
        if single_objective
            [solution, fval, exitflag, output, LB,UB] = ...
                optimizer(objective_function, pc_params, seq);
            % assign dummies for the Pareto fronts
            Paretos = [];

        % call multi-objective optimizer
        elseif multi_objective
            [solution, fval, Paretos, Paretos_fvals, exitflag, output, LB,UB] =...
                optimizer(objective_function, pc_params, seq);
        end

        % call post-processor
        postprocessor_data = [];
        if settings.postprocessing.check && ~cancel_button_pressed() % only when requested
            progress_bar(0, 'Starting post-processor...');
            % run selected post-processor
            index = settings.postprocessing.post_processor;
            [postprocessor_data, solution, fval] = ...
                environment.plugin_info.postprocessors(index).function_handle(...
                @(X) objective_function(X, pc_params), solution, fval, LB, UB);
        end

        %% process results

        % stop the timer
        elapsed_time = toc(start_time);

        % update progress bar
        progress_bar(1, 'Optimization terminated.')

        % evaluate all solutions one more time to obtain all
        % data associated with the solution
        if ~isempty(Paretos) % multi-objective
            % evaluate all Paretos
            Paretos = permute(Paretos, [3,2,1]);
            Paretos_cell = mat2cell(Paretos, ones(size(Paretos,1),1), size(Paretos,2));
            results = cellfun(@(X)patched_conics(seq, X, pc_params), ...
                Paretos_cell, 'uniformoutput', false);
            % gather additional data (if any) from third objective
            if objectives.which(3)
                third_objective_data = cell(size(results));
                for ii = 1:numel(results)
                    [ig,ig, third_objective_data{ii}] = ...
                        feval(objectives.other, MainWin, results{ii});%#ok
                end
            end
            % evaluate the "best" solution among the Paretos
            best_result = find(all( bsxfun(@eq, fval, Paretos_fvals), 2));
            % insert everything in calculation results
            calculation.results.first_order.best.solution             = results{best_result};
            calculation.results.first_order.best.function_value       = fval;
            calculation.results.first_order.Paretos.solutions         = results;
            calculation.results.first_order.Paretos.function_values   = Paretos_fvals;
            if objectives.which(3)
                calculation.results.first_order.best.third_objective_data = ...
                    third_objective_data{best_result};
                calculation.results.first_order.Paretos.third_objective_data = ...
                third_objective_data;
            end
            % also set the types of optimization
            calculation.results.first_order.Paretos.types = settings.optimize.objectives;
            % check feasibility
            violated = cellfun( @(x) x.is_violated, results);
            feasible = ~all(violated);
            calculation.results.first_order.Paretos.infeasible = nnz( violated);
            calculation.results.first_order.Paretos.feasible   = nnz(~violated);

        else % single-objective
            % get all the data
            best_result = patched_conics(seq, solution, pc_params);
            % insert in calculation results
            calculation.results.first_order.best.solution       = best_result;
            calculation.results.first_order.best.function_value = fval;
            calculation.results.first_order.Paretos             = [];
            % check feasibility of the solution
            feasible = ~best_result.is_violated;

        end

        % also save the results from the post-processor
        calculation.results.first_order.best.postprocessor_data = postprocessor_data;

        % get the statistics
        [ephemerides, lambert, lambert_failed, cbf, cbf_failed, rp_failed] = ...
            PC_objective_function(objectives, solution, 'get_stats');

        % insert these data into first-order result
        calculation.results.first_order.elapsed_time     = elapsed_time;
        calculation.results.first_order.lambert_failed   = lambert_failed;
        calculation.results.first_order.lambert_problems = lambert;
        calculation.results.first_order.central_body_problems = cbf;
        calculation.results.first_order.central_body_failed   = cbf_failed;
        calculation.results.first_order.rp_failed        = rp_failed;
        calculation.results.first_order.ephemerides      = ephemerides;
        calculation.results.first_order.algorithm.output = output;
        calculation.results.first_order.algorithm.flag   = exitflag;

        % also save current settings -- these are to be used in the
        % second-order optimziation. Saving the settings here allows the
        % user to modify some settings, while still being able to run the
        % local optimization afterwards and plot the correct results
        calculation.results.first_order.settings = settings;

        % also save the objective function and original bounds -- this is
        % needed for post processors used in combination with multi-objective
        % optimizations
        calculation.results.first_order.objective_function = objective_function;
        calculation.results.first_order.LB = LB;
        calculation.results.first_order.UB = UB;

        % output in case of a normal exit
        varargout{1} = calculation;  % output results
        varargout{2} = true;         % sucessful?
        varargout{3} = false;        % optimize again?

        % if this is a non-batch optimization
        if (nargin > 0) && strcmpi(batch_or_normal, 'normal')
            % Maximum function evaluations or iterations has been reached
            if (exitflag == -1 || exitflag == -2)
                %??? TODO

            % The sequence was REALLY weird - all solutions turned out to
            % be INF or NAN
            elseif (exitflag == -3) && ~cancel_button_pressed()
                % ask whether to optimize again
                button_pressed = questdlg({...
                    'No solution was found to be possible even in theory.';
                    'What would you like to do?'}, 'All solutions are impossible', ...
                    'Try again', 'Give up', 'Give up');
                % optimize again
                if strcmpi(button_pressed, 'Try again')
                    varargout{1} = [];
                    varargout{2} = false;
                    varargout{3} = true;
                    return;
                    % otherwise, just return empty result
                else
                    varargout{1} = [];
                    varargout{2} = false;
                    varargout{3} = false;
                    return;
                end

            % no feasible results
            elseif ~feasible && ~cancel_button_pressed()
                % ask whether to optimize again
                button_pressed = questdlg({'No feasible results were found.';
                    'What would you like to do?'}, 'No results found', ...
                    'Try again', 'Give up', 'Use results anyway', 'Give up');
                % optimize again
                switch lower(button_pressed)
                    case 'try again'
                        varargout{1} = [];
                        varargout{2} = false;   % sucessful?
                        varargout{3} = true;    % optimize again?
                    case 'use results anyway'
                        % output is the default one
                    case 'give up'
                        varargout{1} = [];
                        varargout{2} = false;   % sucessful?
                        varargout{3} = false;   % optimize again?
                end

            % normal exit
            else
                progress_bar(1, 'All done.'), pause(0.25)
            end
        end % if batch or normal

        % We're done!
        % (although not required, give a return  here for clarity)
        return;

    end % first order

    %% Second-order optimization

    if second_order

         %% optimization & post-processing

         % initialize
        progress_bar(0, 'Optimizing...')                 % update progress bar
        pause(0.25), start_time = tic;                   % start timer
        PC_objective_function([],[], 'reset_stats');     % reset persistent variables
        objectives.which = [1 0 0];                      % use only maximum mass objective

        % redefine objective function
        objective_function = @(X, params) PC_objective_function(...
            objectives, patched_conics(seq, X, params));

        % call optimizer
        [solution, fval, exitflag, output] = ...
            optimizer(objective_function, pc_params, seq);

        % call post-processor
        postprocessor_data = [];
        if settings.postprocessing.check && ~cancel_button_pressed() % only when requested
            progress_bar(0, 'Starting post-processor...');
            % TODO: continue here
            [postprocessor_data, ...
             solution, ...
             fval] = post_processor(@(X) objective_function(X, pc_params),...
                                    solution,...
                                    fval,...
                                    LB, UB);
        end

        %% process output

        % stop the timer
        elapsed_time = toc(start_time);

        % update progress bar
        progress_bar(1, 'Optimization terminated.')

        % evaluate solution one more time to obtain all data
        best_result = patched_conics(seq, solution, pc_params);
        % insert in calculation results
        calculation.results.second_order.best.solution       = best_result;
        calculation.results.second_order.best.function_value = fval;
        % check feasibility of the solution
        feasible = ~best_result.is_violated;

        % also save the results from the post-processor
        calculation.results.second_order.best.postprocessor_data = postprocessor_data;

        % get the statistics
        [ephemerides, lambert, lambert_failed, ~, ~, rp_failed] = ...
            PC_objective_function(objectives, solution, 'get_stats');

        % insert these data into first-order result
        calculation.results.second_order.elapsed_time     = elapsed_time;
        calculation.results.second_order.lambert_failed   = lambert_failed;
        calculation.results.second_order.lambert_problems = lambert;
        calculation.results.second_order.rp_failed        = rp_failed;
        calculation.results.second_order.ephemerides      = ephemerides;
        calculation.results.second_order.algorithm.output = output;
        calculation.results.second_order.algorithm.flag   = exitflag;

        % also save the function handle to the solver, so we don't have to
        % re-deinfe it in GENERATE_OUTPUT()
        calculation.results.second_order.solver = objective_function;

        % output in case of a normal exit
        varargout{1} = calculation;  % output results
        varargout{2} = true;         % sucessful?
        varargout{3} = false;        % optimize again?

        % Maximum function evaluations or iterations has been reached
        if (exitflag == -1 || exitflag == -2)
            %??? TODO

        % The sequence was REALLY weird - all solutions turned out to
        % be INF or NAN
        elseif (exitflag == -3) && ~cancel_button_pressed()
            % ask whether to optimize again
            button_pressed = questdlg({'No solution was found to be possible even in theory.';
                'What would you like to do?'}, 'All solutions are impossible', ...
                'Try again', 'Give up', 'Give up');
            % optimize again
            if strcmpi(button_pressed, 'Try again')
                varargout{1} = [];
                varargout{2} = false;
                varargout{3} = true;
                return;
                % otherwise, just return empty result
            else
                varargout{1} = [];
                varargout{2} = false;
                varargout{3} = false;
                return;
            end

            % no feasible results
        elseif ~feasible && ~cancel_button_pressed()
            % ask whether to optimize again
            button_pressed = questdlg({'No feasible results were found.';
                'What would you like to do?'}, 'No results found', ...
                'Try again', 'Give up', 'Use results anyway', 'Give up');
            % optimize again
            switch lower(button_pressed)
                case 'try again'
                    varargout{1} = [];
                    varargout{2} = false;   % sucessful?
                    varargout{3} = true;    % optimize again?
                case 'use results anyway'
                    % output is the default one
                case 'give up'
                    varargout{1} = [];
                    varargout{2} = false;   % sucessful?
                    varargout{3} = false;   % optimize again?
            end

            % normal exit
        else
            progress_bar(1, 'All done.'), pause(0.25)
        end


    end % second order

    % We're done! Yaaay!

end % MGA

%% MGA -- helper functions

% objective/constraint wrapper function for first-order optimization
function varargout = PC_objective_function(objectives, result, stats)

    % initialize persistent variables
    persistent EPHEMERIDES LAMBERT LAMBERT_FAILED RP_FAILED   % statistics
    persistent CENTRAL_BODY CENTRAL_BODY_FAILED               % more statistics

    % initialize
    costs = []; constraints = [];

    % assign initial values (needed for first call) or reset
    if isempty(EPHEMERIDES) || ((nargin == 3) && strcmpi(stats, 'reset_stats'))
        EPHEMERIDES  = 0;   LAMBERT_FAILED = 0;
        LAMBERT      = 0;   RP_FAILED      = 0;
        CENTRAL_BODY = 0;   CENTRAL_BODY_FAILED = 0;
        if ((nargin == 3) && strcmpi(stats, 'reset_stats')), return, end
    end

    % return the statistics
    if (nargin == 3) && strcmpi(stats, 'get_stats')
        varargout{1} =  EPHEMERIDES;     varargout{4} =  CENTRAL_BODY;
        varargout{2} =  LAMBERT;         varargout{5} =  CENTRAL_BODY_FAILED;
        varargout{3} =  LAMBERT_FAILED;  varargout{6} =  RP_FAILED;
        return
    end

    % keep track of the statistics
    LAMBERT        = LAMBERT + result.lambert.success + result.lambert.impossible;
    LAMBERT_FAILED = LAMBERT_FAILED + result.lambert.failure;
    RP_FAILED      = RP_FAILED + result.rp.failure;
    EPHEMERIDES    = EPHEMERIDES + result.ephemerides;
    CENTRAL_BODY   = CENTRAL_BODY + result.central_body.success + result.central_body.impossible;
    CENTRAL_BODY_FAILED = CENTRAL_BODY_FAILED + result.central_body.failure;


    %% minimum Delta_V / maximum drymass

    if objectives.which(1)
        % cost is the negative endmass
        cost = -result.endmass;
        % append all constraints
        rv = result.violations;
        constraint = [...
            rv.C3_violation
            rv.total_DV_violation
            rv.GAM_DV_violation(:)
            rv.rp_violation(:)
            rv.GAM_Vinf_violation(:)
            rv.Vinf_target_violation
            rv.mindist_violation(:)
            rv.TOF_violation
            rv.CBF_violation(:)
            rv.miss_distance(:)
            rv.other_constraints];
        % append costs & constraints
        costs       = [costs, cost];
        constraints = [constraints; constraint];
    end

    %% minimum TOF

    if objectives.which(2)
        % total travel time
        % NOTE: exclude the initial time
        cost = sum(result.tfs(:, 2:end), 2);
        % constraints
        constraint = result.violations.TOF_violation;
        % append costs & constraints
        costs       = [costs, cost];
        constraints = [constraints; constraint];
    end

    %% custom

    if objectives.which(3)
        % calculate the cost & constraints user-defined function
        [cost, constraint] = feval(objectives.other, MainWin, result);
        % append costs & constraints
        costs       = [costs, cost];
        constraints = [constraints; constraint];
    end

    %  output
    varargout{1} = costs;       % objective function value
    varargout{2} = constraints; % all constraints are [c], inequality constraints
    varargout{3} = [];          % [ceq] = [];
    varargout{4} = result;      % (used for post-processing etc.)

end % patched conics objective/constraint function

% small wrapper function for ephemerides generators
function statevec = ephemerides_generation(body, time, model, environment, settings)

    % initialize persistents
    persistent STATES SETTINGS STATIC HAVE_DE405_MEX HAVE_DE421_MEX

    % initialize them on first call
    if (nargin == 5)
         STATES   = model.states;
         SETTINGS = settings;
         STATIC   = model.static;
         
         HAVE_DE405_MEX = environment.have_DE405_MEX;
         HAVE_DE421_MEX = environment.have_DE421_MEX;
         
         return;
    end

    % initialize
    statevec = zeros(1, 6);
    
    % first check static bodies
    if STATIC{1}(body)
        statevec = STATIC{2}{body}; 
        return;
    end

    switch find(SETTINGS.optimize.ephemerides(1),1,'first')
        
        % JPL/DE405
        case 1
            if HAVE_DE405_MEX

                % JPL has different body order; convert
                if body==1, JPL_body = 11;
                else JPL_body = body-1; end

                % JPL expects JD instead of days past J2000
                time_corrected = time + 2451544.5;

                statevec = JPL_DE405_ephemerides(time_corrected, JPL_body);
            else
                statevec = JPL_DE405_ephemerides(time, STATES(body, :));
            end
         
        % JPL/DE421
        case 2
            if HAVE_DE421_MEX

                % JPL has different body order; convert
                if body==1, JPL_body = 11;
                else JPL_body = body-1; end

                % JPL expects JD instead of days past J2000
                time_corrected = time + 2451544.5;

                statevec = JPL_DE421_ephemerides(time_corrected, JPL_body);
                
            else
                error([mfilename ':GUI_inconsistency'], [...
                      'JPL421 selected while MEX file not available; ',...
                      'this indicates a bug in the GUI.']);                
            end
          
        % Kepler (direct)
        case 3
            statevec = two_body_ephemerides(body, time, 0, model,constants);
            % ??? TODO
            
        % quintic best-fit    
        case 4
            statevec = quintic_best_fit(time, body, model,constants);
            % ??? TODO
    end
    
end % ephemerides generation

%% MGA -- main optimization

% global/local optimizer and post-processor
function varargout = optimizer(obj, obj_params, seq)

    %% Initialize

    % get MainWin handle
    global MainWin

    % get relevant data
    settings    = getappdata(MainWin, 'settings'   );
    calculation = getappdata(MainWin, 'calculation');

    % amount of swingby's to perform
    num_swingbys = numel(seq)-2;
    % total amount of targets visited
    num_targets = num_swingbys + 1;

    % handy definitions
    high_thrust  = settings.propulsion.selected(1);
    first_order  = isempty(calculation.results.first_order.best.solution);
    second_order = ~first_order;
    single_objective = nnz(...
        [settings.optimize.objectives.max_mass
        settings.optimize.objectives.min_tof
        settings.optimize.objectives.other.use]) == 1;
    multi_objective  = ~single_objective;

    %% Global optimization
    if first_order

        %% set bounds

        % =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        % FORMAT (Without DSMs):
        %
        %    [t0,...
        %     TOF(1)  longway(1)  k2(1)   m(1)   leftbranch(1),...
        %     TOF(2)  longway(2)  k2(2)   m(2)   leftbranch(2),...
        %     ...
        %     TOF(N)  longway(N)  k2(N)   m(N)   leftbranch(N)]
        %
        % with TOF > 0 the times of flight, longway = +- 1 the
        % corresponding long-way solution, and m = +-1 * max_m
        % the maximum amount of full revolutions to try
        % =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        % initialize LB/UB
        LB = zeros(1, 5*num_targets+1);  UB = LB;

        % include initial/final time
        LB(1) = settings.launch_window_days_past_J2000(1);
        UB(1) = settings.launch_window_days_past_J2000(2);

        % include times of flight, as selected by the user
        TOF_LB_settings = zeros(num_swingbys+1,1);
        TOF_UB_settings = TOF_LB_settings;
        for ii = 1:num_swingbys
            TOF_LB_settings(ii) = settings.GAM.TOF_LB(settings.GAM.body(ii),ii); end
        for ii = 1:num_swingbys
            TOF_UB_settings(ii) = settings.GAM.TOF_UB(settings.GAM.body(ii),ii); end
        TOF_LB_settings(end) = settings.target.TOF_LB(settings.target.body);
        TOF_UB_settings(end) = settings.target.TOF_UB(settings.target.body);
        LB(2:5:end) = TOF_LB_settings;      UB(2:5:end) = TOF_UB_settings;

        % allow long way solutions
        % +1 = short way   -1 = long way
        % setting: 1: short  2: long  3: optimize
        shortlong_settings = ...
            [settings.GAM.short_or_longway(1:num_swingbys); settings.target.short_or_longway];
        lb_longway = ones(num_targets, 1);        % default: short
        ub_longway = ones(num_targets, 1);
        lb_longway(shortlong_settings == 2) = -1; % some are long
        ub_longway(shortlong_settings == 2) = -1;
        lb_longway(shortlong_settings == 3) = -1; % optimize
        % but NOT the upper bound in that case ^_^
        LB(3:5:end) = lb_longway;    UB(3:5:end) = ub_longway;

        % k2-limits
        LB(4:5:end) = 0.01;  UB(4:5:end) = 1; % <- is this OK? I couldn't think of any other limit...

        % complete revolutions [m]
        revolutions = [settings.GAM.complete_revolutions - 1;
            settings.target.complete_revolutions - 1];
        revolutions = revolutions(1:num_targets);
        revolutions_ub = revolutions;  revolutions_lb = revolutions;
        revolutions_ub(revolutions == 3) = 3; % optimize

        revolutions_lb(revolutions == 3) = 0;
        LB(5:5:end) = revolutions_lb;
        UB(5:5:end) = revolutions_ub;
        % +1 = take left-branch   -1 = take right-branch
        if any(UB(5:5:end) ~= 0)
            LB(6:5:end) = -1;   UB(6:5:end) = 1;
        else
            LB(6:5:end) = 0;    UB(6:5:end) = 0;
        end

        %% optimize globally

        % adjust scope
        obj_params.scope = 'global';

        % population size:
        % 100 individuals per independent variable for high-thrust propulsion,
        % 75 for low-thrust propulsions because of the higher computational
        % demands for the low-thrust Lambert targeters
        popsize = (100 - ~high_thrust*25) * numel(LB(LB~=UB));

        % Nothing might have to be done (rare cases of course)
        if all(LB == UB)
            % progress bar
            progress_bar(1, 'Equal bounds - nothing to do.');
            % assign all output
            fval         = obj(LB, obj_params);   varargout{nargout-0} = UB;
            varargout{1} = LB;                    varargout{nargout-1} = LB;
            varargout{2} = fval;
            if multi_objective
                varargout{3} = repmat({LB}, [1, 1, popsize]);
                varargout{4} = repmat(fval, popsize, 1);
                varargout{5} = 1;
                varargout{6}.message = 'Lower and upper bounds were equal - nothing to do.';
            elseif single_objective
                varargout{3} = 1;
                varargout{4}.message = 'Lower and upper bounds were equal - nothing to do.';
            end
            % and return
            return
        end

        % use GODLIKE optimizer
        if (multi_objective || settings.optimize.global.optimizer(1))
            % get current options
            options = settings.optimize.global.optimizer_settings{1};
            % ??? TODO: Better option here
            if isempty(options.NumStreams)
                options.NumStreams = 2; end
            % use only ONE stream for multi-objective optimizations
            if multi_objective
                options.NumStreams = 1; end
            % Don't use PSO for multi-objective optimizations; the current
            % implementation doesn't work so well...
            if multi_objective
                options.algorithms = {'GA';'DE'}; end
            % implement "automatic" popsize
            if isempty(options.popsize)
                options.popsize = popsize; end
            % max. number of function evaluations (needed in the output
            % function)
            maxFunEvals = options.MaxFunEvals;
            % some non-settable options
            options.TolIters = 15;
            options.ReinitRatio = 0.01*single_objective;
            options.AchieveFunVal = -settings.launch.payload_mass;
            options.QuitWhenAchieved =  true;
            options.OutputFcn = @outputFcn;
            options.ConstraintsInObjectiveFunction = 2;
            % and optimize
            [solution, funval, varargout{3:nargout-2}] = ...
                GODLIKE(@(X)obj(X, obj_params), LB, UB, [], options);

        % use MINIMIZE() (FMINSEARCH = NELDERMEAD)
        elseif settings.optimize.global.optimizer(2)
            % get current options
            options = settings.optimize.global.optimizer_settings{2};
            % implement "automatic" popsize
            if isempty(options.popsize)
                options.popsize = popsize; end
            % max. number of function evaluations (needed in the output
            % function)
            maxFunEvals = options.MaxFunEvals;
            % insert non-settable options
            options.OutputFcn = @outputFcn;
            options.ConstraintsInObjectiveFunction = 2;
            options.Algorithm = 'fminsearch';
            % and optimize
            [solution, funval, varargout{3:nargout-2}] = ...
                minimize(@(X)obj(X, obj_params), [], [],[], [],[], LB,UB, [], options);

        % use MINIMIZE() (FMINLBFGS)
        elseif settings.optimize.global.optimizer(3)
            % get current options
            options = settings.optimize.global.optimizer_settings{3};
            % implement "automatic" popsize
            if isempty(options.popsize)
                options.popsize = popsize; end
            % max. number of function evaluations (needed in the output
            % function)
            maxFunEvals = options.MaxFunEvals;
            % insert non-settable options
            options.OutputFcn = @outputFcn;
            options.ConstraintsInObjectiveFunction = 2;
            options.Algorithm = 'fminlbfgs';
            % and optimize
            [solution, funval, varargout{3:nargout-2}] = ...
                minimize(@(X)obj(X, obj_params), [], [],[], [],[], LB,UB, [], options);

        end

        %% optimize locally (after convergence of the global routine)

        % only run local optimizer for single-objective optimizations,
        % and when the cancel button's not been pressed
        if single_objective && ~cancel_button_pressed()
            % change scope
            obj_params.scope = 'local';
            % Check to see if we still have function evaluations left
            maxFunEvals = max(maxFunEvals - varargout{4}.funcCount, 2500);
            % if not, exit IF. otherwise, optimize locally
            if (maxFunEvals > 0)

                progress_bar(0.9, 'Running local optimizer...')

                % set options
                options = optimset(...
                    'display'    , 'off',...        % make sure it's OFF
                    'algorithm'  , 'active-set',... % faster, but less robust (active-set gives NaN/inf often)
                    'FunValCheck', 'on',...         % so, check this!
                    'MaxFunEvals', maxFunEvals,...  % or [maxFunEvals] function evaluations
                    'OutputFcn'  , @outputFcn);     % show intermediate masses in progress bar

                options.ConstraintsInObjectiveFunction = 2;

                % FMINCON
                if settings.optimize.local.optimizer(1)
                    % active-set gives NaN/inf quite often
                    try
                        [solution, funval] = fmincon(@(X)obj(X, obj_params), ...
                            solution,[],[],[],[],LB,UB,@constraints_wrapper,options);
                    % in those cases, use interior-point
                    catch%#ok
                        try % this one too might fail; BARRIER() complains about inf-constraint values
                            options = optimset(options,...
                                'FunValCheck', 'off',...
                                'algorithm'  , 'interior-point');
                            [solution, funval] = fmincon(@(X)obj(X, obj_params), ...
                                solution,[],[],[],[],LB,UB,@constraints_wrapper,options);
                        catch%#ok...so use MINIMIZE(); FMINCON() has utterly failed us
                            options.ConstraintsInObjectiveFunction = 2;
                            options.Algorithm = 'fminsearch';
                            [solution, funval] = minimize(@(X)obj(X, obj_params),...
                                solution, [],[], [],[], LB,UB, [], options);
                        end
                    end

                % Constrained Nelder-Mead
                elseif settings.optimize.local.optimizer(2)
                    options.Algorithm = 'fminsearch';
                    [solution, funval] = minimize(@(X)obj(X, obj_params),...
                        solution, [],[], [],[], LB,UB, [], options);

                % Constrained quasi-Newton (L)BFGS
                elseif settings.optimize.local.optimizer(3)
                    options.Algorithm = 'fminlbfgs';
                    [solution, funval] = minimize(@(X)obj(X, obj_params), ...
                        solution, [],[], [],[], LB,UB, [], options);
                end
            end

        % otherwise, just say the optimization's finished
        else
            progress_bar(1, 'Optimization completed.')
        end

        % assign output arguments
        varargout{1} = solution;
        varargout{2} = funval;
        varargout{nargout-1} = LB; % used in post-processing
        varargout{nargout-0} = UB; % used in post-processing

    end % if first order

    %% High-accuracy local optimization

    if second_order
% TODO : not at all well tested ...

        %% initialize

        % copy ORIGINAL settings
        settings = calculation.results.first_order.settings;

        % redefine high-thrust etc.
        high_thrust = settings.propulsion.selected(1);

        % change scope
        obj_params.scope = 'local';

        % parameters for the optimization
maxFunEvals = 1e4;


        %% define second-order initial estimate and new bounds

        % best result from first-order optimziation
        best_result = calculation.results.first_order.best.solution;
        % initialize initial estimate
        initial_estimate = [best_result.t0; best_result.tfs(:)];

initial_estimate(2) = 8;
initial_estimate(3) = initial_estimate(3) - 8;

        % high-thrust
        if high_thrust
            % insert all departure velocity vectors
            V_init = best_result.V_departure.';
%             initial_estimate = [initial_estimate; V_init(:)];

initial_estimate = [initial_estimate; best_result.DeltaVs(1); V_init(:)];


        % low-thrust
        else
            % Sims & Flanagan
            if settings.optimize.local.method(1)
                % ??? TODO
            % collocation
            elseif settings.optimize.local.method(2)
                % ??? TODO
            end
        end

        % extract original upper and lower bounds
        LB_1st = calculation.results.first_order.LB;
        UB_1st = calculation.results.first_order.UB;

        % set new bounds
        LB = -inf(size(initial_estimate));    UB = +inf(size(initial_estimate)); % no real limits
        LB(1) = LB_1st(1);                    UB(1) = UB_1st(1);                 % launch date
        tofLB = LB_1st(2:5:end);              tofUB = UB_1st(2:5:end);
        LB((0:length(tofUB)-1)+2) = tofLB;    UB((0:length(tofUB)-1)+2) = tofUB; % times of flight

        %% second-order optimization

        % adjust progress bar
        progress_bar(0, 'Running local optimizer...')

        % set options
        options = optimset('display'    , 'off',...         % make sure it's OFF
                           'algorithm'  , 'active-set',...  % faster, but less robust (active-set gives NaN/inf often)
                           'FunValCheck', 'on',...          % so, check this!
                           'MaxFunEvals', maxFunEvals,...   % or [maxFunEvals] function evaluations
                           'OutputFcn'  , @outputFcn);      % show intermediate masses in progress bar
        options.ConstraintsInObjectiveFunction = 2;

        % FMINCON
        if settings.optimize.local.optimizer(1)
            % active-set gives NaN/inf quite often
            try
                [solution, funval] = fmincon(@(X)obj(X, obj_params), ...
                    initial_estimate,[],[],[],[],LB,UB,@constraints_wrapper,options);
                % in those cases, use interior-point
            catch%#ok
                try % this one too might fail; BARRIER() complains about inf-constraint values
                    options = optimset(options,...
                        'FunValCheck', 'off',...
                        'algorithm'  , 'interior-point');
                    [solution, funval] = fmincon(@(X)obj(X, obj_params), ...
                        initial_estimate,[],[],[],[],LB,UB,@constraints_wrapper,options);
                catch%#ok...so use MINIMIZE(); FMINCON() has utterly failed us
                    options.ConstraintsInObjectiveFunction = 2;
                    options.Algorithm = 'fminsearch';
                    [solution, funval] = minimize(@(X)obj(X, obj_params),...
                        initial_estimate, [],[], [],[], LB,UB, [], options);
                end
            end
        % Constrained Nelder-Mead
        elseif settings.optimize.local.optimizer(2)
            options.Algorithm = 'fminsearch';
            [solution, funval] = minimize(@(X)obj(X, obj_params),...
                initial_estimate, [],[], [],[], LB,UB, [], options);
        % Constrained quasi-Newton (L)BFGS
        elseif settings.optimize.local.optimizer(3)
            options.Algorithm = 'fminlbfgs';
            [solution, funval] = minimize(@(X)obj(X, obj_params), ...
                initial_estimate, [],[], [],[], LB,UB, [], options);
        end

        % assign output arguments
        varargout{1} = solution;
        varargout{2} = funval;
        varargout{nargout-1} = LB; % used in post-processing
        varargout{nargout-0} = UB; % used in post-processing

    end % if first/second order

    %% Helper functions

    % small wrapper function to get the nonlinear constraint values into
    % the FMINCON() optimization routine
    %
    % All the other optimizers I wrote understand how to get constraint
    % function values from within objective functions, thus evaluating the
    % objective function only *once* to get both objective function values
    % and constraint function values. But FMINCON() can't do that, so it
    % needs to evaluate PATHCED_CONICS() *twice* to get the constraint
    % values.
    %
    % This is indeed an obvious imperfection with regard to performance,
    % but we don't really care about that here -- FMINCON() on
    % PATHCED_CONICS() is only used for either a few iterations after the
    % global routine, or in one of the post-processors, in which case the
    % largest computational burder (by far) is in the objective function,
    % and definitely not here.
    function [c, ceq] = constraints_wrapper(X)
        [ignored, c, ceq] = obj(X, obj_params); %#ok<ASGLU>
    end

    % This nested function serves as a general output function for
    % all optimizers. It updates the progress bar and keeps track
    % of the optimization statistics. It also causes the
    % optimization to stop when the cancel-button's been pressed, or
    % certain obvious termination conditions apply
    function stop = outputFcn(X, optimValues, state)%#ok

        % only stop if the cancel button has been pressed. But stop
        % IMMEDIATELY, e.g., don't display values from the last iteration
        stop = cancel_button_pressed();
        if stop
            return; end

        % First correct the sleepy idiot who wrote the output function
        % sections of FMINSEARCH() and FMINCON()
        if isfield(optimValues, 'funccount')
            % funCCCCCCCount (CAPITAL!)
            optimValues.funcCount = optimValues.funccount;
            optimValues           = rmfield(optimValues, 'funccount');
        end

        % initialize
        show_progress = false;

        % global optimizers
        if first_order

            % discriminate between the different optimizers
            if any(strcmpi(state, {'iter';'interrupt'}))

                % Constrained FMINSEARCH(), FMINCON() and Constrained FMINLBFGS()
                if (strcmpi(obj_params.scope, 'local') && rand < 1/5) ||...
                   (strcmpi(obj_params.scope, 'global') && ...
                   ~strcmpi(state, 'interrupt') && isfield(optimValues, 'procedure'))
                    % get display values
                    endmass = abs(optimValues.fval);
                    if isfield(optimValues, 'constrviolation')
                        if isfield(optimValues.constrviolation, 'nonlin_ineq')
                            % MINIMIZE()
                            violation = optimValues.constrviolation.nonlin_ineq{2};
                        else
                            % FMINCON()
                            violation = optimValues.constrviolation;
                        end
                    else
                        violation = 0;
                    end
                    % update progress
                    if strcmpi(obj_params.scope, 'local')
                        progress = 0.9 + 0.1*optimValues.funcCount/maxFunEvals;
                    else
                        progress = max([...
                            optimValues.funcCount/maxFunEvals *0.9, ...
                            endmass/settings.launch.launch_mass *0.9]);
                    end
                    show_progress = true;
                    % stop immediately if mass equals given launch mass and
                    % constraint violation is zero
                    if abs(endmass - settings.launch.launch_mass) < 1e-1 &&...
                       violation == 0;
                        stop = true;
                    end
                end

                % GODLIKE; single-objective optimizations
                if isfield(optimValues, 'type') && ...
                   strcmpi(optimValues.type, 'single-objective')
                    show_progress = true;
                    % modify the display string
                    violation = optimValues.best_fval_constraint_violation;
                    endmass   = abs(optimValues.best_fval);
                    % update progress
                    if (violation == 0)
                        progress = max([optimValues.funcCount/maxFunEvals *0.9
                                        endmass/settings.launch.launch_mass *0.9]);
                    elseif ~isfinite(endmass)
                        endmass  = 0;
                        progress = 0;
                    elseif (violation > 0)
                        %pen_endmass = endmass*(1 - (exp(violation) - 1));
                        progress = max([...
                            optimValues.funcCount/maxFunEvals *0.9
                            endmass/settings.launch.launch_mass *0.9]);
                            %pen_endmass*/settings.launch.launch_mass *0.9
                    else
                        progress = 0.90;
                    end
                    % stop immediately if mass equals given launch mass and
                    % constraint violation is zero
                    if abs(endmass - settings.launch.launch_mass) < 1e-1 &&...
                       violation == 0;
                        stop = true;
                    end
                end

                % update progress and progress bar
                if show_progress
                    display_string = ['Final payload mass: ', ...
                        num2str(endmass,5), ' kg'];
                    if (violation > 0) || (endmass < settings.launch.payload_mass)
                        display_string = [display_string, ' (INFEASIBLE)'];
                    end
                    progress_bar(progress, display_string);
                end

                % GODLIKE; multi-objective optimizations
                if isfield(optimValues,'type') && ...
                        strcmpi(optimValues.type, 'multi-objective')
                    progress = ...
                        optimValues.number_of_nondominated_solutions/optimValues.popsize;
                    display_string = sprintf('%3.2f%% of solutions non-dominated',...
                        100*optimValues.number_of_nondominated_solutions/optimValues.popsize);
                    progress_bar(progress, display_string);
                end % if (multi-objective)

            end % if (iter)
        end % first order (global) output

        % local optimizations
        if second_order
            %??? TODO
        end % second order (local) output

    end % outputFcn

end % optimizer & post-processor
