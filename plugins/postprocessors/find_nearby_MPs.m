function varargout = find_nearby_MPs(varargin)
% FIND_NEARBY_MPS         Find nearby MP's
%
%
%

    %% return info & function handles 
    
    % build the infostructure
    pp_info.name = 'Find mp''s close to the trajectory';
    pp_info.GUI_function_handle = @(value) GUI_interaction(value);
    pp_info.plot_function_handle = @(varargin) plot_fcn(varargin{:});
    pp_info.function_handle = @(obj, x0, fval, LB, UB)...
                              pp_fcn(obj, x0, fval, LB, UB);
    
    % return this info on initial (zero-argument) call
    if (nargin == 0), varargout{1} = pp_info; return; end
    
    %% What should we do when it's selected in the GUI? 
    
    function GUI_interaction(value)
        
        %% initialize
        
        % extract required data
        global MainWin arrival_tab sequence_tab
        settings    = getappdata(MainWin, 'settings'   );
        handles     = getappdata(MainWin, 'handles'    );
        environment = getappdata(MainWin, 'environment');
        
        % make some abbreviations
        at = handles.tab(arrival_tab);
        st = handles.tab(sequence_tab);
        % see if the MPCORB database has been loaded. If not,
        % prompt the user whether to load it or not
        if ~settings.model(4)
            % dialog
            Button_pressed = questdlg(...
                {'This post-processor requires the MPCORB-database to be loaded.'
                ' '
                'Do you wish to enable that option?'
                ' '
                'NOTE: because the MPCORB contains data on ~400.000 minor planets, '
                'this post-processor can easily become be a very costly operation. '
                'Make sure all the mission specific MP-pruning functions (in '
                './MISSION_SPECIFIC/) are as strict as possible for the current '
                'mission.'}, 'MPCORB needs to be loaded', 'Yes','Cancel','Yes');
            % decide what to do
            switch lower(Button_pressed)
                case 'yes'
                    % check the checkbox
                    MP_check_handle = st.selected_model(4);
                    set(MP_check_handle, 'value', 1);
                    % change the setting
                    settings.model(4) = true;
                    setappdata(MainWin, 'settings', settings);
                    % run its callback function
                    callbacks('MPs_check', false, MP_check_handle);
                case 'cancel'
                    % reset the post-processor to previous value;
                    set(at.post_processing.post_processor(2),...
                        'value', settings.postprocessing.post_processor);
                    % reset the appropriate fields
                    if isfield(at.post_processing, 'additional_controls') && ...
                       isfield(at.post_processing.additional_controls, 'panel')
                        set(at.post_processing.additional_controls.panel(:), ...
                            'visible', 'off');
                    end
                    % don't forget to change its setting
                    settings.postprocessing.post_processor = 1;
                    setappdata(MainWin, 'settings', settings);
                    % and return
                    return
            end
        end   
        
        %% draw the GUI components
                
        % insert additional settings
        
        % change setting
        settings.postprocessing.post_processor = value;
        
        % Distance-dependent or constant threshold
        handles.tab(arrival_tab).post_processing.additional_controls.panel(1) = uipanel(...
            'parent'    , at.post_processing_panel,...
            'position'  , [0.01 0.2 0.98 0.6],...
            'visible'   , 'off',...
            'bordertype', 'none'); 
        at = handles.tab(arrival_tab); % it's changed, so redefine it

        % explanatory text
        uicontrol(...
            'parent'  , at.post_processing.additional_controls.panel(1),...
            'units'   , 'normalized',...
            'style'   , 'text',...
            'horizontalalignment', 'left',...
            'position', [0.01 0.77 0.98 0.23],...
            'string'  , {'In order to find minor planets "close" to the optimal trajectory,'
                         'you have to specify what "close" means. '
                         'The routine will accept MP''s that come closer to the trajectory than '
                         'a certain threshold distance:'});

        % threshold selection        
        threshold_group = uibuttongroup(...
            'parent'  , at.post_processing.additional_controls.panel(1),...
            'position', [0 0.2 1 0.6],...
            'bordertype', 'none',...
            'SelectionChangeFcn', @(varargin){callbacks('find_nearby_MPs_threshold', varargin{:})
                modify_settings('change_single_setting', 'postprocessing.threshold.check',...
                'handles.tab(arrival_tab).post_processing.additional_controls.threshold_radio',...
                varargin{:})});
            % constant threshold 
            handles.tab(arrival_tab).post_processing.additional_controls.threshold_radio(1) = uicontrol(...
                'style'   , 'radiobutton', ...
                'parent'  , threshold_group,...
                'units'   , 'normalized',...
                'position', [0.01 0.76 1 0.15],...
                'value'   , settings.postprocessing.threshold.check(1),...
                'string'  , 'use constant threshold:');
            handles.tab(arrival_tab).post_processing.additional_controls.constant_threshold.text(1) = uicontrol(...
                'style'   , 'text',...
                'units'   , 'normalized',...
                'parent'  , threshold_group,...
                'enable'  , 'off',...
                'string'  , 'threshold = ',...
                'position', [0.2 0.5 0.3 0.15]);
            handles.tab(arrival_tab).post_processing.additional_controls.constant_threshold.value = uicontrol(...
                'style'   , 'edit',...
                'Background', environment.colors.edit_bgcolor,...
                'string'  , settings.postprocessing.constant_threshold.distance,...
                'enable'  , 'off',...
                'units'   , 'normalized',...
                'parent'  , threshold_group,...
                'position', [0.5 0.52 0.15 0.15],...
                'callback', @(varargin) modify_settings('change_single_setting',...
                    'postprocessing.constant_threshold.distance', [], varargin{:}));
            handles.tab(arrival_tab).post_processing.additional_controls.constant_threshold.text(2) = uicontrol(...
                'style'   , 'text',...
                'units'   , 'normalized',...
                'enable'  , 'off',...
                'parent'  , threshold_group,...
                'string'  , '(AU)',...
                'position', [0.7 0.5 0.15 0.15]);

            % threshold depends on S/C apocenter
            handles.tab(arrival_tab).post_processing.additional_controls.threshold_radio(2) = uicontrol(...
                'style'   , 'radiobutton', ...
                'parent'  , threshold_group,...
                'units'   , 'normalized',...
                'position', [0.01 0.35 1 0.15],...
                'value'   , settings.postprocessing.threshold.check(2),...
                'string'  , 'threshold depends on S/C apocenter (ra):');                
            handles.tab(arrival_tab).post_processing.additional_controls.variable_threshold.text(1) = uicontrol(...
                'style'   , 'text',...
                'units'   , 'normalized',...
                'parent'  , threshold_group,...
                'string'  , 'threshold = ',...
                'position', [0.0 0.1 0.20 0.15]);
            handles.tab(arrival_tab).post_processing.additional_controls.variable_threshold.value(1) = uicontrol(...
                'style'   , 'edit',...
                'Background', environment.colors.edit_bgcolor,...
                'string'  , settings.postprocessing.variable_threshold.distance(1),...
                'units'   , 'normalized',...
                'parent'  , threshold_group,...
                'position', [0.2 0.12 0.15 0.15],...
                'callback', @(varargin) modify_settings('change_single_setting',...
                    'postprocessing.variable_threshold.distance(1)', [], varargin{:}));
            handles.tab(arrival_tab).post_processing.additional_controls.variable_threshold.text(2) = uicontrol(...
                'style'   , 'text',...
                'units'   , 'normalized',...
                'parent'  , threshold_group,...
                'string'  , '+',...
                'position', [0.36 0.1 0.05 0.15]);
            handles.tab(arrival_tab).post_processing.additional_controls.variable_threshold.value(2) = uicontrol(...
                'style'   , 'edit',...
                'Background', environment.colors.edit_bgcolor,...
                'string'  , settings.postprocessing.variable_threshold.distance(2),...
                'units'   , 'normalized',...
                'parent'  , threshold_group,...
                'position', [0.41 0.12 0.15 0.15],...
                'callback', @(varargin) modify_settings('change_single_setting',...
                    'postprocessing.variable_threshold.distance(2)', [], varargin{:}));
            handles.tab(arrival_tab).post_processing.additional_controls.variable_threshold.text(3) = uicontrol(...
                'style'   , 'text',...
                'units'   , 'normalized',...
                'parent'  , threshold_group,...
                'string'  , '* ra /',...
                'position', [0.55 0.1 0.15 0.15]);
            handles.tab(arrival_tab).post_processing.additional_controls.variable_threshold.value(3) = uicontrol(...
                'style'   , 'edit',...
                'Background', environment.colors.edit_bgcolor,...
                'string'  , settings.postprocessing.variable_threshold.distance(3),...
                'units'   , 'normalized',...
                'parent'  , threshold_group,...
                'position', [0.70 0.12 0.15 0.15],...
                'callback', @(varargin) modify_settings('change_single_setting',...
                    'postprocessing.variable_threshold.distance(3)', [], varargin{:}));
            handles.tab(arrival_tab).post_processing.additional_controls.variable_threshold.text(4) = uicontrol(...
                'style'   , 'text',...
                'units'   , 'normalized',...
                'parent'  , threshold_group,...
                'string'  , '(AU)',...
                'position', [0.84 0.1 0.15 0.15]);

            % also try to minimize distances
            uicontrol(...
                'style'   , 'check',...
                'units'   , 'normalized',...
                'parent'  , at.post_processing.additional_controls.panel(1),...
                'position', [0.01, 0.1, 1, 0.1],...
                'string'  , 'Also try to minimize the distances to the nearby MP''s',...
                'value'   , settings.postprocessing.optimize_distance.check,...
                'callback', @(varargin)...
                    modify_settings('change_single_setting', ...
                    'postprocessing.optimize_distance.check', [], varargin{:}));
        
        % show/hide additional GUI-components
% TODO: this ain't generalized yet...
set(at.post_processing.additional_controls.panel(2:end), 'visible', 'off');
set(at.post_processing.additional_controls.panel(1), 'visible', 'on');
        
        % save the handles
        setappdata(MainWin, 'handles', handles);        
    end % GUI function  
    
    %% What should we do when any results are to be plotted? 
    
    function plot_fcn(varargin)
        
        % get appropriate data
        global MainWin output_tab trajectory_tab post_processing
        calculation = getappdata(MainWin, 'calculation');    
        model       = getappdata(MainWin, 'model'      );
        handles     = getappdata(MainWin, 'handles'    );
        constants   = getappdata(MainWin, 'constants'  );
                
        % rename some constants for brevity
        AU  = constants.AU;
        muC = constants.mu_central;
        R = calculation.results.first_order.best.postprocessor_data;
        ot = handles.tab(output_tab);
        trajectory_axes = ot.tab(trajectory_tab).trajectory_pane(1);        
        postprocessor_infopane = ot.tab(post_processing).postprocessor_infopane;
        postprocessor_title    = ot.tab(post_processing).postprocessor_title;
        
        % make sure we're plotting in the trajectory plot
        set(MainWin, 'currentaxes', trajectory_axes);
    
% TODO: centralize this!
textcolor = 'w';
txtoffset = 0.1;
        
        % set table title
        set(postprocessor_title, 'string', [num2str(R.nearby_MPs), ...
            ' Minor planet''s found close to the optimal trajectory']);
        % only plot when at least one MP was found
        if (R.nearby_MPs > 0)
            % logical index to reachable MP's
            reachables = any(R.reachable, 2);
            % extract MP-names and initial times
            names = model.MPs.MP_names(reachables);
            t0s   = model.MPs.epoch_days_past_J2000(reachables);
            % get positions of the reachable MP's
            xs = model.MPs.x0s(reachables,:);
            for ii = 1:size(xs,1)
                % extract encounter times
                times = R.encounter_times(ii,:);
                times = times(isfinite(times));
                % find the associated positions
                states = progress_orbit(times(:)-t0s(ii), xs(ii, :), muC);
                % plot MP immediately
                % ??? TODO: make a nice warp 
                plot3(states(:,1)/AU, states(:,2)/AU, states(:,3)/AU, ...
                    'color', textcolor, ...
                    'marker', '.');
                % also plot its name
                text(states(:, 1)/AU + txtoffset,...
                    states(:, 2)/AU + txtoffset,...% also works when the same MP is found
                    states(:, 3)/AU, names{ii}, ...% multiple times
                    'color'   , textcolor,...
                    'clipping', 'on');
                
            end % nearby-MP plot loop
            
% for now..."it just works"
% ??? TODO: multiples
R.encounter_times(~isfinite(R.encounter_times)) = 0;
R.encounter_times = sum(R.encounter_times,2);
R.min_dist(~isfinite(R.min_dist)) = 0; 
R.min_dist = sum(R.min_dist,2);
R.rel_speed(~isfinite(R.rel_speed)) = 0; 
R.rel_speed = sum(R.rel_speed,2);
            
            % convert & collect data
            data = {};
            for ii = (1:R.nearby_MPs)
                encounter_date = create_date(R.encounter_times(ii));
                data = [data; {names{ii}, R.rel_speed(ii), R.min_dist(ii), ...
                    R.min_dist(ii)/constants.AU, '', encounter_date(1:14)}];%#ok
            end
            
            % define table
            % min. dist, rel. velocity, ...
            set(postprocessor_infopane, ...
                'columnformat', {'char', 'numeric', 'numeric', 'numeric', [], 'char'},...
                'ColumnName'  , {'Designation/Name', 'Relative speed |(km/s)', ...
                'closest approach | (km)', 'closest approach | (AU)',...
                '', 'Encounter date'}, ...
                'data'        , data,...
                'columnwidth' , {145, 105, 105, 105, 50, 105});
            
            % if not MP's were found, clear the table
        else
            set(postprocessor_infopane, 'data', {});
        end
        
        % create date-string from given days past J2000.0
        % TODO: redundant copy-paste from generate output...
        function date = create_date(days)
            monthstring = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug',...
                'Sep' 'Oct' 'Nov' 'Dec'};
            [y, M, d, h, m] = days2date(days);
            day    = append_suffix(d); 
            if numel(day) == 3, day = [' ',day]; end
            year   = num2str(y);
            hour   = num2str(h);  if numel(hour) == 1, hour = ['0',hour]; end
            minute = num2str(m);  if numel(minute) == 1, minute = ['0', minute]; end
            date   = [monthstring{M}, ' ', day, ', ', year, ', ', hour, ':', ...
                minute, ' (GMT)'];            
        end % create date
        
        % create proper counter suffix
        % TODO: redundant copy-paste from generate output...
        function out_string = append_suffix(integer)
            % create counter string
            inp_string = num2str(integer);
            if any(strcmp(inp_string,{'11';'12';'13'}))
                out_string = [inp_string,'th'];
            else
                switch inp_string(end)
                    case '1', out_string  = [inp_string,'st'];
                    case '2', out_string  = [inp_string,'nd'];
                    case '3', out_string  = [inp_string,'rd'];
                    otherwise, out_string = [inp_string,'th'];
                end
            end
        end % append suffix
        
    end % plot function
    
    %% What is the actual post-processing routine?
    
    %if (settings.postprocessing.post_processor == 2)
    function [postprocessor_data, solution, fval] = ...
            pp_fcn(obj, x0, fval, LB, UB)
        
        %% initialize
        
        global MainWin 
        settings  = getappdata(MainWin, 'settings' );        
        constants = getappdata(MainWin, 'constants');
        model     = getappdata(MainWin, 'model'    );
        
        % initially, the new solution is the given initial value
        postprocessor_data = [];
        solution = x0;
        
        % high thrust or low thrust?
        high_thrust = settings.propulsion.selected(1);
        low_thrust  = ~high_thrust;
        
        %% build initial close-proximity MP-database
        
        % progress bar
        progress_bar(0, 'Finding MP''s close to the optimal trajectory...')
        % (initial) options for the minimum distance routines
        if settings.postprocessing.threshold.check(1)
            mindist_options = struct('threshold', ...
                settings.postprocessing.constant_threshold.distance*constants.AU);
        elseif settings.postprocessing.threshold.check(2)
            mindist_options = struct('threshold', ...
                settings.postprocessing.variable_threshold.distance*constants.AU);
        end
        
        % find nearby MPs
        postprocessor_data = find_nearby_MPs(x0);
        
        % update progress bar
        if (postprocessor_data.nearby_MPs == 0)
            progress_bar(1, 'No nearby MP''s were found.');
        elseif postprocessor_data.nearby_MPs == 1
            progress_bar(1, [num2str(postprocessor_data.nearby_MPs),...
                ' MP found within threshold distance.']);
        else
            progress_bar(1, [num2str(postprocessor_data.nearby_MPs),...
                ' MP''s found within threshold distance.']);
        end
        
        %% minimize distances (if requested and at least 1 nearby MP was found)
        
        if settings.postprocessing.optimize_distance.check && ...
          (postprocessor_data.nearby_MPs > 0)
            
            %% initialize
            
            % find indices to nearby MP's (only once)
            indices = find(any(postprocessor_data.reachable, 2));
            % initial minimum distances in AU
            % ??? TODO: multiples! This formulation "just works"!
            D0 = postprocessor_data.min_dist; D0(~isfinite(D0)) = 0; D0 = sum(D0,2);
            % progress bar
            progress_bar(0, 'Attempting to minimize all distances...')
            % maximum function evaluations for all optimizers
            maxfuneval = 500;
            % adjust bounds; allow 10 day variation around the travel times
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
            LB(1) = x0(1) - 5;                UB(1) = x0(1) + 5;
            LB(2:5:end) = x0(2:5:end) - 5;    UB(2:5:end) = x0(2:5:end) + 5;
            
            %% run local optimizer
            
            % FMINCON
            if settings.optimize.local.optimizer(1)
                % set options
                options = optimset(...
                    'maxfunevals', maxfuneval,...
                    'display'    , 'off');
                % define proper objective & constraint fucntions
                objective  = @(X) mindist_to_MPs_objective_function(X, 'cost');
                constraint = @(X) mindist_to_MPs_objective_function(X, 'constraint');
                % and optimize
                solution = ...
                    fmincon(objective, x0, [],[],[],[], LB,UB, constraint, options);
            % Constrained Nelder-Mead
            elseif settings.optimize.local.optimizer(2)
                % set options
                options = optimset(...
                    'MaxFunEvals', maxfuneval,...
                    'display'    , 'off');
                options.ConstraintsInObjectiveFunction = 2;
                % define proper objective & constraint fucntions
                objective  = @(X) mindist_to_MPs_objective_function(X, 'both');
                % and optimize
                solution = ...
                    optimize(objective, x0, [],[],[],[], LB,UB, [],[], options, 'fminsearch');
            % Constrained quasi-Newton (L)BFGS
            elseif settings.optimize.local.optimizer(3)
                % set options
                options = optimset(...
                    'MaxFunEvals', maxfuneval,...
                    'display'    , 'off');
                options.ConstraintsInObjectiveFunction = 2;
                % define proper objective & constraint fucntions
                objective  = @(X) mindist_to_MPs_objective_function(X, 'both');
                % and optimize
                solution = ...
                    optimize(objective, x0, [],[],[],[], LB,UB, [],[], options, 'fminlbfgs');
            end
            
            %% evaluate new solution
            
            [fval, postprocessor_data] = ...
                mindist_to_MPs_objective_function(solution, 'final');
                        
        end % if 
        
        %% Helper functions

        % find nearby MPs - either build initial database, or use existing
        % database to compare distances with new trial solution
        function [postprocessor_data, cost, eq_constraints] = find_nearby_MPs(X, indices)

            % evaluate patched conics
            [cost, dummy, eq_constraints, result] = obj(X);%#ok
            % convert transfer and initial times to something more convenient
            times = cumsum([result.t0, result.tfs]);

            % select mode of operation:
            % build initial database, or do an iteration in the optimizer
            init = false; if (nargin == 1), init = true; indices = 1:numel(model.MPs.as); end

            % initialize loop
            reachable       = false(size(model.MPs.as(indices),1),numel(result.tfs));% logical
            encounter_times = zeros(size(reachable));% double
            min_dist        = encounter_times;       % double
            rel_speed       = encounter_times;       % double

            % loop through all trajectories in the result
            for jj = 1:numel(result.tfs)

                % build initial database. this is a lengthy process, so we
                % need an output function
                if init
                    mindist_options.outputFcn = @(state,values) mindist_outputFcn(...
                        (jj-1)/numel(result.tfs), jj/numel(result.tfs),...
                        ['Leg ', num2str(jj), ': '], state,values);
                    % otherwise, optimize the distance. No output function needed
                else
                    mindist_options.outputFcn = [];
                end

                % put all MP's at the proper initial positions
                % NOTE: don't forget to convert times to seconds!
                t0 = times(jj); tend = times(jj+1);
                MPs_new_Ms = model.MPs.M0s(indices) + ...
                    model.MPs.ns(indices).*(t0-model.MPs.epoch_days_past_J2000(indices))*86400;
                MPs_new_theta0s = eM2theta(model.MPs.es(indices), MPs_new_Ms);

                % run the minimum distance routine
                % (high-thrust post-processor differs from low-thrust post-processor)
                if high_thrust
                    % convert all trajectories therein to orbital elements
                    sc_elements = cart2kep([result.r_departure, result.V_departure], ...
                        constants.mu_central, 'theta');
                    % run routine
                    [reachable(:,jj),encounter_times(:,jj),min_dist(:,jj),rel_speed(:,jj)] = ...
                        minimum_distance_conics(sc_elements(jj,:), ...
                        [model.MPs.as(indices), model.MPs.es(indices), model.MPs.is(indices), ...
                        model.MPs.Omegas(indices), model.MPs.omegas(indices), MPs_new_theta0s], ...
                        t0, tend, constants.mu_central, mindist_options);
                elseif low_thrust
                    % run routine
                    [reachable(:,jj),encounter_times(:,jj),min_dist(:,jj),rel_speed(:,jj)] = ...
                        minimum_distance_exposins(result.exposins(jj,:), result.R{jj},...
                        [model.MPs.as(indices), model.MPs.es(indices), model.MPs.is(indices), ...
                        model.MPs.Omegas(indices), model.MPs.omegas(indices), MPs_new_theta0s], ...
                        t0, tend, constants.mu_central, mindist_options);
                end
            end % minimum-distance for loop

            % gather data
            reachables = any(reachable,2);
            postprocessor_data.reachable        = reachable;
            postprocessor_data.min_dist         = min_dist(reachables,:);
            postprocessor_data.rel_speed        = rel_speed(reachables,:);
            postprocessor_data.encounter_times  = encounter_times(reachables,:);
            postprocessor_data.nearby_MPs       = nnz(reachables);

        end

        % objective / constraint function for the minimization of the overall
        % distances to the nearby MP's
        function varargout = mindist_to_MPs_objective_function(X, mode)

            % initialize
            Dnew = inf(postprocessor_data.nearby_MPs,1);

            % evaluate patched conics & minimum distance routine,
            % this time with the indices of only the nearby MP's
            [new_postprocessor_data, cost, eq_constraints] = find_nearby_MPs(X, indices);

            % rename those that are reachable
            reachable = (sum(new_postprocessor_data.reachable, 2) > 0);

            % new minimal distances
            % ??? TODO: multiples! This formulation "just works"!
            Dnew_temp = new_postprocessor_data.min_dist;
            Dnew_temp(~isfinite(Dnew_temp)) = 0; Dnew_temp = sum(Dnew_temp,2);

            % assign values to proper indices
            Dnew(reachable) = Dnew_temp; % the rest's just [inf]

            % inequality constraints: (Dnew - D0 <= 0) to ensure all distances
            % are less than the initial distances found
            ineq_constraints = Dnew - D0;

            % new cost function: original cost (max. mass ONLY, no
            % multi-objective here) PLUS the sum of the distances in AU
            costs = cost(1)*1000 - sum(Dnew(isfinite(Dnew)));

            % assign output arguments
            switch lower(mode)
                case 'cost' % cost function only
                    varargout{1} = costs;            % -mass - sum(distances)
                case 'constraint' % constraint function only
                    varargout{1} = ineq_constraints; % (c   <= 0)
                    varargout{2} = eq_constraints;   % (ceq == 0)
                case 'both' % both cost and constraint function
                    varargout{1} = costs;            % -mass - sum(distances)
                    varargout{2} = ineq_constraints; % (c   <= 0)
                    varargout{3} = eq_constraints;   % (ceq == 0)
                case 'final' % last evaluation to get all the proper values
                    varargout{1} = cost;             % new function value at the new solution
                    % modify new postprocessor data
                    new_reachable = false(size(model.MPs.as,1),...
                        size(new_postprocessor_data.reachable,2));
                    new_reachable(indices, :) = new_postprocessor_data.reachable;
                    new_postprocessor_data.reachable = new_reachable;
                    % and assign to output
                    varargout{2} = new_postprocessor_data; % new postprocessor data
            end

        end % minimum distance to nearby MP's cost/constraint function

        % output function for minimum-distance-to-MP's-routine
        % A relatively primitive way to show the progress of the routine that
        % computes the minimum distances to nearby MP's (post-processor)
        function stop = mindist_outputFcn(min_progress, max_progress, string, state, values)
            % evaluate state ofcancel button
            if (nargout == 1), stop = cancel_button_pressed; end
            % progress indicator range
            progress_range = max_progress - min_progress;
            % apses check
            if strcmpi(state, 'apsescheck')
                progress_bar(min_progress, [string, 'apses check: ', ...
                    num2str(values.accepted),' of ',num2str(values.candidates), ...
                    ' MP''s accepted']);
            end
            % overlap check
            if strcmpi(state, 'overlapcheck')
                progress_bar(min_progress + progress_range*values.iteration/values.candidates, ...
                    [string, 'overlap check: ', ...
                    num2str(values.accepted),' of ',num2str(values.candidates), ' MP''s accepted']);
            end
        end % mindist outputFcn
        
    end % "find nearby MP's" post-processor
    
end % find nearby MP's
