% generate all output & graphics
function generate_output(embedding)

    %% Initialize
    % ==========================================================================

    % get all globals
    global MainWin output_tab
    global Pareto_tab trajectory_tab central_body_speed post_processing
    global BATCH_optimization optimization_statistics

    % get application data
    calculation = getappdata(MainWin, 'calculation');
    constants   = getappdata(MainWin, 'constants'  );
    model       = getappdata(MainWin, 'model'      );
    handles     = getappdata(MainWin, 'handles'    );
    environment = getappdata(MainWin, 'environment');

    % some abbreviations
    ot = handles.tab(output_tab);

    % use the settings from the optimization, NOT the current settings
    settings = calculation.results.first_order.settings;

    % enable the output tab-button
    set(handles.tabbutton(output_tab), 'enable', 'on');

    % rename for clarity & brevity
    firstorder  = calculation.results.first_order;
    secondorder = calculation.results.second_order;
    highthrust  = settings.propulsion.selected(1);
    ionengine   = settings.propulsion.selected(2);
    Solarsail   = settings.propulsion.selected(3);
    singleobjective = isempty(firstorder.Paretos);
    multiobjective  = ~singleobjective;

    % rename results
    results = firstorder;

    % some renaming and initializations
    AU  = constants.AU;
    muC = constants.mu_central;
    seq = results.best.solution.seq(:).';
    monthstring = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
    status = '';
    previous_point_trajectory = [];
    previous_point_Paretos    = [];

% get options
% TODO: FUTURE WORK
np = 500;
transferorbitcolor   = [.2 .2 .2];
transfersectioncolor = 'g';
transfersectioncolor_bad = 'r';
hyperbolafac         = .95;
GAMbody_orbitcolor = [.4, .4, .8];
txtoffset        = 0.1;
settings.plot_options.resolution = 'low';

    % defaults to all embedded
    Pareto_axes            = ot.tab(Pareto_tab).Pareto_pane(1);
    Pareto_infopane        = ot.tab(Pareto_tab).Pareto_pane(2);
    trajectory_axes        = ot.tab(trajectory_tab).trajectory_pane(1);
    trajectory_infopane    = ot.tab(trajectory_tab).trajectory_pane(2);
    central_speed_axes     = ot.tab(central_body_speed).central_body_speed_pane;
    postprocessor_infopane = ot.tab(post_processing).postprocessor_infopane;
    postprocessor_title    = ot.tab(post_processing).postprocessor_title;

    MainFig(1:3)           = MainWin;

    trajectory_external    = false;
    Paretos_external       = false;
    central_external       = false;
    BATCH_external         = false;

    % first list optimization statistics
    print_optimization_statistics;

    % re-define when embedding is different
    if strcmpi(embedding, 'embedded')
        % do nothing

    elseif strcmpi(embedding, 'separate_trajectory')
        trajectory_external  = true;
        MainFig(2)           = figure('renderer', 'opengl');
        trajectory_axes      = axes('parent', MainFig(2));
        trajectory_infopane  = [];

    elseif strcmpi(embedding, 'separate_Paretofront')
        Paretos_external = true;
        MainFig(1)       = figure('renderer', 'opengl');
        Pareto_axes      = axes('parent', MainFig(1));
        Pareto_infopane  = [];
    end

    % call the appropriate nested function
    % second order plots
    if ~isempty(secondorder.best.solution)
        % second-order, high-thrust
        if highthrust
            high_second;
        % second-order, low-thrust
        elseif (ionengine || Solarsail)
            low_second;
        end

    % first-order plots
    else
        % first order, single-objective, high-thrust
        if singleobjective && highthrust
            high_first_single(false);

        % first order, multi-objective, high-thrust
        elseif multiobjective && highthrust
            high_first_multi;

        % first order, single-objective, low-thrust
        elseif singleobjective && (ionengine || Solarsail)
            low_first_single(false);

        % first order, multi-objective, low-thrust
        elseif multiobjective && (ionengine || Solarsail)
            low_first_multi;
        end
    end

    %% Plot single-objective low-thrust trajectories
    % ==========================================================================

    function low_first_single(done)

        % Initialization & infopane
        % ======================================================================

        % don't re-plot the trajectory if the Paretos are external
        if Paretos_external, return, end

        % by default, we haven't been here before
        if (nargin == 0), done = false; end

        % enable trajectory & central velocity button
        set(ot.tab(trajectory_tab).button, 'enable', 'on')
        set(ot.tab(central_body_speed).button, 'enable', 'on')

        % first put some global information about the solution in the infopane
        % (for embedded plots)
        if ~trajectory_external
            solution_info = create_solution_info;
            set(trajectory_infopane, 'string', solution_info);
            % also set the UserData field (for the "reset" option)
            set(trajectory_infopane, 'UserData', solution_info);
        end

        % Plot trajectory
        % ======================================================================

        % first plot all the GAM-bodies and their orbits
        plot_GAMbody_trajectories(done);

        % rename for clarity
        R = results.best.solution;

        % plot all exposins in current solutionvector
        for j = 1:numel(seq)-1

            % extract parameters from this exposin
            [k0, k1, k2, phi, dth] = deal(R.exposins(j,1), R.exposins(j,2), ...
                R.exposins(j,3), R.exposins(j,4), R.exposins(j,7));

            % define path
            r  = @(th) k0.*exp(k1.*sin(k2.*abs(th) + phi));
            th = (linspace(0, dth, np)).';

            % apply coordinate transformation
            [x,y] = pol2cart(th, r(abs(th))); z = zeros(size(x));
            xyz = [x,y,z]*R.R{j}/AU;

            % determine whether the GAM was feasible. If it was, plot
            % green, otherwise, plot red
            good_bad_color = transfersectioncolor;
            if (R.violations.total_DV_violation > 0 ||...
                R.violations.TOF_violation > 0)
                % all swingbys are invalid
                good_bad_color = transfersectioncolor_bad;
            else
                % launch
                if (j == 1)
                    if (R.violations.C3_violation > 0)
                        good_bad_color = transfersectioncolor_bad; end
                % swingbys
                elseif (R.violations.GAM_DV_violation(j-1) > 0 ||...
                        R.violations.rp_violation(j-1) > 0)
                    good_bad_color = transfersectioncolor_bad;
                end
            end

            % and plot
            plot3(xyz(:,1), xyz(:,2), xyz(:,3), good_bad_color, 'linewidth', 2);
        end

        % correct the possibly distorted axes
        if ~trajectory_external, view(0, -90);
        else view(120, 30);
        end
        axis equal

        % Plot speed w.r.t. central body
        % ======================================================================

        % don't plot the central speed if an external plot of the
        % trajectory is requested
        if trajectory_external
            return; end

        % set the current axes to the central speed axes
        set(MainFig(3), 'currentaxes', central_speed_axes);
        cla reset
        grid on
        hold on

        % define times
        times = cumsum([results.best.solution.t0, results.best.solution.tfs]);

        % helper functions
        r     = @(th) k0*exp(k1*sin(k2*th+phi));
        thdot = @(th) sqrt( muC./r(th).^3 .* ...
            (1 + k1^2*k2^2*cos(k2*th+phi).^2 + k1*k2^2*sin(k2*th+phi)) );
        rdot  = @(th) thdot(th).*k1*k2*cos(k2*th+phi).*r(th);

        % plot
%??? TODO

        % set axis labels, title etc.
        ylabel('speed [km/s]')
        title(['Spacecraft speed w.r.t. ', model.CentralBody{1}])
        % reset limits and tick marks
        xticks = get(central_speed_axes, 'xtick');
        set(central_speed_axes, 'xtick', linspace(xticks(1), xticks(end), 10));
        xticks = get(central_speed_axes, 'xtick');
        yticks = get(central_speed_axes, 'ytick');
        for ii = 1:numel(xticks)
            date{ii} = create_date(xticks(ii), 'short'); end%#ok
        set(central_speed_axes, 'xticklabel', []);
        text(xticks-.5*(xticks(2)-xticks(1)), ...
            repmat(yticks(1)-.45*(yticks(2)-yticks(1)),length(xticks),1),...
            date, 'rotation', 20);

    end % low first single

    %% Plot multi-objective low-thrust trajectories
    % ==========================================================================

    function low_first_multi
        % first plot the "most efficient" point
        low_first_single();
        % don't plot Paretos if the trajectory is externalized
        if trajectory_external
            return; end
        % plot the Pareto's
        plot_pareto_front();
    end % low first multi


    %% Plot single-objective high-thrust trajectories
    % ==========================================================================

    % Plots trajectories from a first-order, single-objective,
    % high-thrust solution to the MGA-problem
    function high_first_single(done)

        % Initialization & infopane
        % ======================================================================

        % don't re-plot the trajectory if the Paretos are external
        if Paretos_external
            return; end

        % by default, we haven't been here before
        if (nargin == 0)
            done = false; end

        % enable trajectory & central velocity button
        set(ot.tab(trajectory_tab).button, 'enable', 'on')
        set(ot.tab(central_body_speed).button, 'enable', 'on')
        if ~done
            set(ot.tab(Pareto_tab).button, 'enable', 'off'), end

        % first put some global information about the solution in the infopane
        % (for embedded plots)
        if ~trajectory_external
            solution_info = create_solution_info;
            set(trajectory_infopane, 'string', solution_info);
            % also set the UserData field (for the "reset" option)
            set(trajectory_infopane, 'UserData', solution_info);
        end

        % Plot trajectory
        % ==========================================================================

        % first plot all the GAM-bodies and their orbits
        plot_GAMbody_trajectories(done);

        % rename some variables
        r_departure = results.best.solution.r_departure;
        V_departure = results.best.solution.V_departure;
        r_arrival   = results.best.solution.r_arrival;
        V_arrival   = results.best.solution.V_arrival;
        m           = results.best.solution.m;

        % compute orbital elements of the transfer trajectories
        orbitals1 = cart2kep([r_departure, V_departure], muC, 'theta');
        orbitals2 = cart2kep([r_arrival, V_arrival], muC, 'theta');
        orbitals  = [orbitals1, orbitals2(:, end)];

        % plot all trajectories
        for j = 1:size(orbitals, 1)

            % get orbital data for the current transfer
            a = orbitals(j, 1);  O = orbitals(j, 4);
            e = orbitals(j, 2);  o = orbitals(j, 5);
            i = orbitals(j, 3);

            % define thetas and positions
            thstart = orbitals(j, 6);
            thend   = orbitals(j, 7);
            Vdep    = V_departure(j, :);

            % treat elliptic orbits and hyperbolic trajectories seperately
            ell = (e < 1);   hyp = (e > 1);

            % elliptic orbit
            if ell
                % put thetas in 0 <= theta <= 2pi
                thstart = mod(thstart, 2*pi);     thend = mod(thend, 2*pi);
                % determine direction of motion
                xtrans1 = kep2cart(a, e, i, O, o, thstart + 1e-6, muC, 'theta');
                xtrans2 = kep2cart(a, e, i, O, o, thstart - 1e-6, muC, 'theta');
                % define transfer orbits and transfer sections
                smallint = abs(thend - thstart);
                bigint   = 2*pi - smallint;
                % this method to decide in which direction to plot is a
                % "nightly build"...I'm sure it can be done simpler :)
                if m(j) == 0
                    if (sign(xtrans1(1) - xtrans2(1)) == sign(Vdep(1)))
                        if (thstart > thend)
                            thtransfer    = linspace(0, bigint, np)' + thstart;
                            thnontransfer = linspace(0, smallint, np)' + thend;
                        else
                            thtransfer    = linspace(0, smallint, np)' + thstart;
                            thnontransfer = linspace(0, bigint, np)' + thend;
                        end
                    else
                        if (thstart > thend)
                            thnontransfer = linspace(0, bigint, np)' + thstart;
                            thtransfer    = linspace(0, smallint, np)' + thend;
                        else
                            thnontransfer = linspace(0, smallint , np)' + thstart;
                            thtransfer    = linspace(0, bigint, np)' + thend;
                        end
                    end
                else
                    thtransfer    = linspace(0, 2*pi, np)';
                    thnontransfer = zeros(size(thtransfer));
                end

            % hyperbolic trajectory
            elseif hyp
                smallint = abs(thend - thstart);
                if (thend < thstart)
                    thtransfer     = linspace(0, smallint, np)' + thend;
                    thnontransfer1 = linspace(-hyperbolafac*acos(-1/e), thend, np)';
                    thnontransfer2 = linspace(thstart, hyperbolafac*acos(-1/e), np)';
                else
                    thtransfer     = linspace(0, smallint, np)' + thstart;
                    thnontransfer1 = linspace(-hyperbolafac*acos(-1/e), thstart, np)';
                    thnontransfer2 = linspace(thend, hyperbolafac*acos(-1/e), np)';
                end
            end

            % clone orbitals for vectorization
            elements = repmat([a, e, i, O, o], np, 1);

            % full transfer orbit/trajectory (non-transfer part)
            if ell
                xyznontrans1 = kep2cart([elements, thnontransfer], muC, 'theta');
                xyznontrans1 = xyznontrans1(:, 1:3)/AU;
            elseif hyp
                xyznontrans1 = kep2cart([elements, thnontransfer1], muC, 'theta');
                xyznontrans2 = kep2cart([elements, thnontransfer2], muC, 'theta');
                xyznontrans1 = xyznontrans1(:, 1:3)/AU;
                xyznontrans2 = xyznontrans2(:, 1:3)/AU;
            end

            % transfer part
            xyztrans = kep2cart([elements, thtransfer], muC, 'theta');
            xyztrans = xyztrans(:, 1:3)/AU;

            % determine whether the GAM was feasible. If it was, plot
            % green, otherwise, plot red
            good_bad_color = transfersectioncolor;
            if (results.best.solution.violations.total_DV_violation > 0 ||...
                    results.best.solution.violations.TOF_violation > 0)
                % all swingbys are invalid
                good_bad_color = transfersectioncolor_bad;
            else
                % launch
                if (j == 1)
                    if (results.best.solution.violations.C3_violation > 0)
                        good_bad_color = transfersectioncolor_bad;
                    end
                % swingbys
                elseif (results.best.solution.violations.GAM_DV_violation(j-1) > 0 ||...
                        results.best.solution.violations.rp_violation(j-1) > 0)
                    good_bad_color = transfersectioncolor_bad;
                end
            end

            % make plots
            nontransfer = plot3(xyznontrans1(:, 1), xyznontrans1(:, 2), xyznontrans1(:, 3));
            if hyp
                nontransfer = [nontransfer; ...
                    plot3(xyznontrans2(:, 1), xyznontrans2(:, 2), xyznontrans2(:, 3))];%#ok
            end
            transfer = plot3(xyztrans(:, 1), xyztrans(:, 2), xyztrans(:, 3));
            set(nontransfer, 'Color', transferorbitcolor);
            set(transfer, 'Color', good_bad_color, 'LineWidth', 2);
        end % for

        % finish trajectory plot
        if ~trajectory_external, view(0, -90);
        else view(120, 30);
        end
        axis equal

        % Plot speed w.r.t. central body
        % ======================================================================

        % don't plot the central speed if an external plot of the
        % trajectory is requested
        if trajectory_external
            return; end

        % set the current axes to the central speed axes
        set(MainFig(3), 'currentaxes', central_speed_axes);
        cla reset, grid on, hold on

        % define times
        times = cumsum([results.best.solution.t0, results.best.solution.tfs]);

        % calculate the central speed as a function of time
        for j = 1:size(orbitals, 1)

            % replicate orbital data for the current transfer
            a  = repmat(orbitals(j, 1), np, 1);
            e  = repmat(orbitals(j, 2), np, 1);
            i  = repmat(orbitals(j, 3), np, 1);
            O  = repmat(orbitals(j, 4), np, 1);
            o  = repmat(orbitals(j, 5), np, 1);

            % time range
            time_range = linspace(times(j), times(j+1), np);

            % convert these times to [theta]'s
            th_range = aet2theta(a, e, time_range(:), times(j), muC);

            % shift to appropriate theta (impossible to automate
            % because coordinate transformation is indifferent to time)
            th_range = th_range - th_range(1) + orbitals(j, 6);

            % get velocities
            [xdot,xdot,xdot, xdot,ydot,zdot] = ...
                kep2cart(a, e, i, O, o, th_range, muC, 'theta');%#ok

            % compute speeds
            speeds = sqrt(sum([xdot, ydot, zdot].^2,2));

            % append everything to plot
            plot(time_range, speeds, 'r', 'linewidth', 2);

            %
            if (j == 1)
                % indicator line
                line([time_range(1), time_range(1)], ...
                    [min(speeds), max(speeds)], 'color', 'k')
                % accompanying text
                text(time_range(1)-5, (max(speeds)+min(speeds))/2, ...
                    ['Launch from ', model.names{seq(1)}], ...
                    'rotation', 90,...
                    'fontweight', 'bold',...
                    'horizontalalignment', 'center',...
                    'clipping', 'on');
            end

            if j ~= size(orbitals, 1)
                % indicator line
                line([time_range(end), time_range(end)], ...
                    [min(speeds), max(speeds)], 'color', 'k')
                Vendprev = speeds(end);
                % accompanying text
                text(time_range(end)-5, (max(speeds)+min(speeds))/2, ...
                    [model.names{seq(j+1)}, ' encounter'], ...
                    'rotation', 90,...
                    'fontweight', 'bold',...
                    'horizontalalignment', 'center',...
                    'clipping', 'on');
            else
                % indicator line
                line([time_range(end), time_range(end)], ...
                    [min(speeds), max(speeds)], 'color', 'k')
                % accompanying text
                text(time_range(end)-5, (max(speeds)+min(speeds))/2, ...
                    ['Arrival at ', model.names{seq(end)}], ...
                    'rotation', 90,...
                    'fontweight', 'bold',...
                    'horizontalalignment', 'center',...
                    'clipping', 'on');
            end

            % set correct length of previous indicator line
            if (j > 1)
                line([time_range(1), time_range(1)], ...
                    [Vendprev, speeds(1)], 'color', 'k')
            end

        end

        % set axis labels, title etc.
        ylabel('speed [km/s]')
        title(['Spacecraft speed w.r.t. ', model.CentralBody{1}])

        % reset limits and tick marks
        xticks = get(central_speed_axes, 'xtick');
        set(central_speed_axes, 'xtick', linspace(xticks(1), xticks(end), 10));
        xticks = get(central_speed_axes, 'xtick');
        yticks = get(central_speed_axes, 'ytick');
        for ii = 1:numel(xticks)
            date{ii} = create_date(xticks(ii), 'short');%#ok
        end
        set(central_speed_axes, 'xticklabel', []);
        text(xticks-.5*(xticks(2)-xticks(1)), ...
            repmat(yticks(1)-.45*(yticks(2)-yticks(1)),length(xticks),1),date, 'rotation', 20)

    end % high first single

    %% Plot multi-objective high-thrust trajectories
    % ==========================================================================

    function high_first_multi
        % first plot the "most efficient" point
        high_first_single(false);
        % don't plot Paretos if the trajectory is externalized
        if trajectory_external, return, end
        % plot the Pareto's
        plot_pareto_front;
    end

    %% Plot Pareto fronts
    % ==========================================================================

    function plot_pareto_front()

        % Get the current x and y limits on that figure
        if ~Paretos_external
            xlim_original_trajectory = get(trajectory_axes, 'xlim');
            ylim_original_trajectory = get(trajectory_axes, 'ylim');
        end

        % enable tab-button
        set(ot.tab(Pareto_tab).button, 'enable', 'on')

        % print statistics & info
        Pareto_info = {...
            ' '
            ['  Feasible Pareto points  : ', num2str(results.Paretos.feasible)]
            ' ';
            ['  Infeasible Pareto points: ', num2str(results.Paretos.infeasible)]
            ' ';
            ''};
        set(Pareto_infopane, 'string', Pareto_info);

        % rename for brevity
        bestfval = results.best.function_value;
        Paretos  = results.Paretos.solutions;
        funvals  = results.Paretos.function_values;
        types    = results.Paretos.types;

        % first objective is ALWAYS the (negative) mass
        funvals(:, 1) = abs(funvals(:, 1));
        bestfval      = abs(bestfval);

        % plot only possible for 2 or 3-objective optimizations
        if size(funvals,2) <= 3

            % initialize axes
            set(0, 'currentfigure', MainFig(1));
            set(MainFig(1), 'currentaxes', Pareto_axes), cla reset, hold on
            hold on

            % endmass is ALWAYS included in the objectives
            labels{1} = 'endmass [kg]';

            % the second objective is time of flight
            if types.min_tof
                labels{2} = 'time of flight [days]';
                % scale to years if all TOF's are larger than than 2 years
                if all(funvals(:, 2) > 365.25*2)
                    labels{2} = 'time of flight [years]';
                    funvals(:, 2) = funvals(:, 2)/365.25;
                    results.Paretos.function_values(:,2) = ...
                        results.Paretos.function_values(:,2)/365.25;
                end
            end

            % the second/third/fourth... objectives are "other" objectives
            if types.other.use
                labels{end+1} = settings.optimize.objectives.other.axis_label; end

            % find those Paretos with constraint violations
            non_violated = cellfun(@(x) ~x.is_violated, Paretos);
            violated     = ~non_violated;
            % don't use logicals (more complicated later on)
            non_violated = find(non_violated);
            violated     = find(violated);

            % to use the mouseover effects, it's easier to plot each
            % point individually and save a handle for each point

            % 2-D fronts
            if size(funvals,2) == 2

                % plot feasible frontmembers
                good_Pareto_points = zeros(numel(non_violated),1);
                for ii = 1:numel(non_violated)
                    good_Pareto_points(ii) = ...
                        plot(funvals(non_violated(ii), 1), ...
                             funvals(non_violated(ii), 2), '.', ...
                             'color', [.9, .2, .9],...'color', [.2, .9, .0],...
                             'MarkerSize', 15);
                end
                % also plot INfeasible frontmembers
                bad_Pareto_points = zeros(numel(violated),1);
                for ii = 1:numel(violated)
                    bad_Pareto_points(ii) = ...
                        plot(funvals(violated(ii), 1), ...
                        funvals(violated(ii), 2), 'rx');
                end
                % emphasize currently selected point
                plot(bestfval(1), bestfval(2), 'k.', 'MarkerSize', 20);
                % set axes labels
                xlabel(labels{1});  ylabel(labels{2});
                % reverse the X-direction
                % (endmass was NEGATIVE during optimizations)
                set(Pareto_axes, 'XDir', 'reverse')
                % scale the axes to focus the display on the feasible points
                if any(non_violated)
                    axis([min(funvals(non_violated, 1)), max(funvals(non_violated, 1)),...
                        min(funvals(non_violated, 2)), max(funvals(non_violated, 2))]);
                % if there are no feasible points, just scale to fit data
                elseif (size(funvals,1) > 1)
                    axis([min(funvals(:, 1)), max(funvals(:, 1)),...
                        min(funvals(:, 2)), max(funvals(:, 2))]);
                end
                % set mouse-click function for all the Pareto points
                set(good_Pareto_points, 'buttondownfcn', @(varargin)...
                    plot_new_trajectory(good_Pareto_points, bad_Pareto_points, varargin{:}));

                % get current limits on axes
                xlim_original = get(Pareto_axes, 'xlim');
                ylim_original = get(Pareto_axes, 'ylim');

                % enable mouseover functions
                % NOTE: these must be an *extention* to the already existing
                % mouseover functions
                if ~Paretos_external
                    set(MainFig(1), ...
                        'WindowScrollWheelFcn', @(varargin) {...
                            scroll_zoom(trajectory_axes, varargin{:})
                            scroll_zoom(Pareto_axes, varargin{:})},...
                        'WindowButtonDownFcn', @(varargin) {...
                            pan_click(trajectory_axes, xlim_original_trajectory, ...
                            ylim_original_trajectory, varargin{:})
                            pan_click(Pareto_axes, xlim_original, ylim_original, varargin{:})},...
                        'WindowButtonUpFcn' , @(varargin) {...
                            pan_release(trajectory_axes, varargin{:})
                            pan_release(Pareto_axes, varargin{:})},...
                        'WindowButtonMotionFcn', @(varargin) {...
                            pan_motion(trajectory_axes, varargin{:})
                            pan_motion(Pareto_axes, varargin{:})
                            Paretos_mouseover(good_Pareto_points, bad_Pareto_points, varargin{:})});
                end

            % 3-D fronts
            else
                %??? TODO
            end

        % for Pareto fronts of more than 3 dimensions, we have
        % to make some sort of smart list, or table
        else
            %??? TODO
        end

        grid on

    end % high first multi


    %% Plot second-order high-thrust trajectory
    % ==========================================================================

    function high_second
        %??? TODO
    end % high second


    %% Plot second-order low-thrust trajectory
    % ==========================================================================

    function low_second
        %??? TODO

        % Was Sims & Flanagan's method or Collocation used?

    end % low second


    %% Helper functions
    % ==========================================================================

    % list optimization statistics
    function print_optimization_statistics

        % enable the tab-button
        set(ot.tab(optimization_statistics).button, 'enable', 'on');

        % statistics for first-order optimization

        % header
        firstorder_header = {...
            ' '; '  Elapsed time';
            ' '; '  Ephemerides calculations';
            ' '; '  Lambert problems solved';
            ' '; '  Lambert problems failed';
            ' '; '  Central body problems solved';
            ' '; '  Central body problems failed';
            ' '; '  Pericenter calculations failed'};

        % handle Pareto's
        if multiobjective
            firstorder_header = [firstorder_header;
                ' '; '  Feasible Pareto points';
                ' '; '  Infeasible Pareto points'];
        end
        firstorder_header = [firstorder_header;
            ' '; '  Algorithm exitflag';
            ' '; '  Algorithm messsage'];
        set(ot.tab(optimization_statistics).firstorder_statistics_pane(1),...
            'string', firstorder_header);

        % separator
        firstorder_separator = repmat(':', size(firstorder_header,1),1);
        firstorder_separator(1:2:end) = ' ';
        set(ot.tab(optimization_statistics).firstorder_statistics_pane(2),...
            'string', firstorder_separator);

        % data
        firstorder_data = {...
            ' '; [num2str(results.elapsed_time), ' sec'];
            ' '; num2str(results.ephemerides);
            ' '; num2str(results.lambert_problems);
            ' '; num2str(results.lambert_failed);
            ' '; num2str(results.central_body_problems);
            ' '; num2str(results.central_body_failed);
            ' '; num2str(results.rp_failed)};
        if multiobjective
            firstorder_data = [firstorder_data;
                ' '; num2str(results.Paretos.feasible);
                ' '; num2str(results.Paretos.infeasible)];
        end
        firstorder_data = [firstorder_data;
            ' '; num2str(results.algorithm.flag);
            ' '; results.algorithm.output.message];
        set(ot.tab(optimization_statistics).firstorder_statistics_pane(3),...
            'string', firstorder_data);


        % statistics for second-order optimization
        % ??? TODO

    end % optimization statistics

    % plot GAM-body trajectories into the trajectory plot
    function plot_GAMbody_trajectories(done)

        % Plot gam-bodies
        % ======================================================================

        % by default, we haven't plotted this one before
        if (nargin == 0)
            done = false; end

        % initialize axes
        set(0, 'currentfigure', MainFig(2));
        set(MainFig(2), 'currentaxes', trajectory_axes), cla reset, hold on
        if trajectory_external
            % for external plots, the trajectory plot already exists.
            % Get the current color of the axes
            current_color = ...
                get(ot.tab(trajectory_tab).trajectory_pane(1), 'color');
            set(trajectory_axes, 'color', current_color)
            textcolor = ~current_color;
        else
            % Default is to use black background for better contrast
            set(trajectory_axes, 'color', 'k')
            textcolor = 'w';
        end

        % extract and rename some parameters
        rps = [results.best.solution.r_departure;
               results.best.solution.r_arrival(end, :)];
        Vps = [results.best.solution.V_departure_GAMbody
               results.best.solution.V_target_GAMbody(end, :)];

        % get the GAM-body orbits
        GAMbody_orbits = cart2kep([rps, Vps], muC, 'theta');

        % build title string for orbitplot
        seqstring = 'Swingby-sequence ';
        bodynames = cell(numel(seq), 1);
        for ii = 1:numel(seq)
            name          = model.names{seq(ii)};
            bodynames(ii) = {name};
            seqstring     = [seqstring, name(1:2)]; %#ok
        end

        % prepare planet orbits
        plth = linspace(0, 2*pi, np);
        a    = repmat(GAMbody_orbits(:, 1), 1, np);
        e    = repmat(GAMbody_orbits(:, 2), 1, np);
        in   = repmat(GAMbody_orbits(:, 3), 1, np);
        o    = repmat(GAMbody_orbits(:, 4), 1, np);
        O    = repmat(GAMbody_orbits(:, 5), 1, np);
        plth = repmat(plth, 1, size(a,1));
        [plx, ply, plz] = kep2cart(a, e, in, O, o, plth, muC, 'theta');
        [plx, ply, plz] = deal(plx/AU, ply/AU, plz/AU);

        % process resolution setting
        if strcmpi(settings.plot_options.resolution, 'extra low')
            step = 8; sphere_size = 10;
        elseif strcmpi(settings.plot_options.resolution, 'low')
            step = 4; sphere_size = 20;
        elseif strcmpi(settings.plot_options.resolution, 'normal')
            step = 2; sphere_size = 50;
        elseif strcmpi(settings.plot_options.resolution, 'high')
            step = 1; sphere_size = 250;
        elseif strcmpi(settings.plot_options.resolution, 'extra high')
            step = 1; sphere_size = 500;
        end

        % create simple sphere to plot GAM-bodies
        [x, y, z] = sphere(sphere_size);

        % initialize handles
        body = zeros(size(plx, 1), 1); % handles to refer to plotted bodies
        body_orbit = body;             % handles to refer to plotted orbits

        % plot every body's orbit
        for j = 1:numel(bodynames)

            % Make exception for the central body
            if strcmpi(bodynames(j), model.CentralBody{1})
                % Print name
                text(rps(j, 1)/AU + txtoffset,...
                     rps(j, 2)/AU + txtoffset,...
                     rps(j, 3)/AU, 'Central body Flyby', ...
                     'color'   , textcolor,...
                     'clipping', 'on');
                 % Draw line from central body to pericenter
                 line([0 rps(j, 1)/AU], [0 rps(j, 2)/AU], [0 rps(j, 3)/AU], ...
                     'color'    , 'y',...
                     'linewidth', 2);
                 % Continue with the next body
                continue;
            end

            % if the current body is a body from the standard model,
            % textures etc. can definitely be used
            if seq(j) <= numel(model.mean_Radii)

                % set body radius
                if j == 6 % ??? TODO: This is only valid for the Sun..
                    Radius = 100 * model.mean_Radii(seq(j))/AU;
                else
                    Radius = 250 * model.mean_Radii(seq(j))/AU;
                end

                % check if the body has a texture
                if size(model.textures, 1) >= seq(j) && ...
                   ~isempty(model.textures{seq(j),1})
                    % retreive departure body's texture
                    texture = model.textures{seq(j), 1};
                    % adjust to resolution setting
                    texture = texture(1:step:end, 1:step:end, :);

                % if the body does not have a texture, just plot
                % a boring gray sphere
                else
                    texture = ones(10,10,3)/2;
                end

                % plot body
                body(j) = warp(...
                    x*Radius + rps(j, 1)/AU, ...
                    y*Radius + rps(j, 2)/AU, ...
                    z*Radius + rps(j, 3)/AU, texture);

                % plot body names
                text(rps(j, 1)/AU + txtoffset,...
                     rps(j, 2)/AU + txtoffset,...
                     rps(j, 3)/AU, bodynames(j), ...
                     'color'   , textcolor,...
                     'clipping', 'on');

                % plot planet orbits
                body_orbit(j) = plot3(plx(j, :), ply(j, :), plz(j, :));
                set(body_orbit(j), 'Color', GAMbody_orbitcolor, 'UserData', 0);

            % otherwise, check if there is a texture
            else
                %??? TODO
                body(j) = NaN;
                body_orbit(j) = NaN;

                % plot body name
                text(rps(j, 1)/AU + txtoffset,...
                    rps(j, 2)/AU + txtoffset,...
                    rps(j, 3)/AU, bodynames(j), ...
                    'color', textcolor,...
                    'clipping', 'on');

            end
        end % for

        % also plot central body
        Radius  = 5 * model.mean_Radii(1)/AU;         % scale radius
        texture = model.textures{1, 1};               % retreive its texture
        texture = texture(1:step:end, 1:step:end, :); % adjust to resolution setting
        warp(x*Radius, y*Radius, z*Radius, texture);  % plot the body

        % include labels & title
        xlabel('X (AU)')
        ylabel('Y (AU)')
        zlabel('Z (AU)');
        title(seqstring);

        % get original axes limits
        xlim_original = get(trajectory_axes, 'xlim');
        ylim_original = get(trajectory_axes, 'ylim');

        % enable zooming with the mouse scroll wheel and panning with
        % clicks (only for embedded figures)
        if ~trajectory_external && ~done
            set(MainFig(2), ...
                'WindowScrollWheelFcn' , @(varargin)...
                    scroll_zoom(trajectory_axes, varargin{:}),...
                'WindowButtonDownFcn'  , @(varargin) ...
                    pan_click(trajectory_axes, xlim_original, ylim_original, varargin{:}),...
                'WindowButtonUpFcn'    , @(varargin)...
                    pan_release(trajectory_axes, varargin{:}),...
                'WindowButtonMotionFcn', @(varargin)...
                    pan_motion(trajectory_axes, varargin{:}));
        end

        %% Plot results from post-processor
        % ======================================================================

        % only if a post-processor was selected
        if settings.postprocessing.check &&...
           ~isempty(results.best.postprocessor_data)
            % first enable button
            set(ot.tab(post_processing).button, 'enable', 'on');
            % then call the appropriate plotting function
            index = settings.postprocessing.post_processor;
            environment.plugin_info.postprocessors(index).plot_function_handle()
        end % post-processor

        %% Plot additional data from third-objective
        % ======================================================================

        % only if a post-processor was selected
        if false && isfield(results.best, 'third_objective_data') &&...
           ~isempty(results.best.third_objective_data)
            % rename for brevity
            R = results.best.third_objective_data;
            % enable button and change its name
            set(ot.tab(post_processing).button, ...
                'enable', 'on',...
                'string', 'Third objective');

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
                    states = progressOrbit(times(:)-t0s(ii), xs(ii, :), muC);
% plot MP immediately
% TODO: ...warp?
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
R.encounter_times(~isfinite(R.encounter_times)) = 0;R.encounter_times = sum(R.encounter_times,2);
R.min_dist(~isfinite(R.min_dist)) = 0; R.min_dist = sum(R.min_dist,2);
R.rel_speed(~isfinite(R.rel_speed)) = 0; R.rel_speed = sum(R.rel_speed,2);

                % convert & collect data
                data = {};
                for ii = (1:R.nearby_MPs)
                    encounter_date = create_date(R.encounter_times(ii), 'full');
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

        end % third objective data

    end % plot_GAMbody_trajectories

    %· · · · · · · · · · · · · ·
    % string building functions
    %· · · · · · · · · · · · · ·

    % create solution information string
    function info_string = create_solution_info

        % default spacing strings
        spaces1   = repmat(' ', 1, 2);
        spaces2   = repmat(' ', 1, 6);
        separator = repmat('-', 1, 35);

        % rename for brevity
        R = results.best.solution;

        % convert the launch date to calendar date
        info_string = {...
            [spaces1, 'Launch from ', model.names{seq(1)}]
            separator
            [spaces2, create_date(R.t0, 'full')]
            [spaces2, 'Launch mass: ', num2str(R.masses(1), 5), ' kg']
            [spaces2, 'C3: ', num2str(R.C3, 5), ' km' char(178) '/s' char(178)]};
        if (R.violations.C3_violation > 0)
            info_string = [info_string; {[spaces1, '(C3 CONSTRAINT VIOLATED)']; ' '}];
        else
            info_string = [info_string; ' '];
        end

        % Calculate the encounter dates in days past J2000.0
        GAM_days = cumsum([results.best.solution.t0, results.best.solution.tfs]);
        GAM_days = GAM_days(2:end);
        % convert them to Calendar dates
        for ii = 1:numel(GAM_days)-1

            info_string = [info_string; {...
                [spaces1, model.names{seq(ii+1)}, ' encounter'];
                separator
                [spaces2, create_date(GAM_days(ii), 'full')]
                [spaces2, 'mass before: ', num2str(R.masses(ii  ), 5), ' kg']
                [spaces2, 'mass after : ', num2str(R.masses(ii+1), 5), ' kg']
                [spaces2, 'turnangle: ', num2str(R.turnangle(ii)*180/pi,5), char(176)]
                [spaces2, 'DV @ pericenter: ', num2str(R.DeltaVs(ii),5), ' km/s']
                [spaces2, 'additional DV: ', num2str(R.additional_DeltaV(ii),5), ' km/s']
                [spaces2, 'Closest approach: ', num2str(R.closest_approach(ii),5), ' km']
                ' '}];%#ok
        end
         info_string = [info_string; {...
                [spaces1, 'Arrival at ', model.names{seq(end)}];
                separator
                [spaces2, create_date(GAM_days(end), 'full')]
                [spaces2, 'arrival V-inf: ', num2str(sqrt(R.Vinf_target*R.Vinf_target.'), 5), ' km/s']
                ' '}];

        % Delta-V
        info_string = [info_string; {...
            [spaces1, 'Total Delta-V: ', num2str(R.DeltaVtot, 5), ' km/s']}];
        % final mass
        info_string = [info_string; {...
            [spaces1, 'Final payload mass: ', num2str(R.masses(end), 5), ' kg']}];
        % travel time
        info_string = [info_string; {...
            [spaces1, 'Total travel time: ', num2str(sum(R.tfs)/365.25, 3), ' years']}];

    end % parse dates

    % create date-string from given days past J2000.0
    function date = create_date(days, len)
        [y, M, d, h, m] = days2date(days);
        day    = append_suffix(d); if numel(day) == 3, day = [' ',day]; end
        year   = num2str(y);
        hour   = num2str(h);  if numel(hour) == 1, hour = ['0',hour]; end
        minute = num2str(m);  if numel(minute) == 1, minute = ['0', minute]; end
        if strcmpi(len, 'full')
            date = [monthstring{M}, ' ', day, ', ', year, ', ', hour, ':', minute, ' (GMT)'];
        elseif strcmpi(len, 'short')
            date = [monthstring{M}, ' ', year];
        end
    end % create date

    % create proper counter suffix
    % (as in 2 -> '2nd',  3 -> '3rd', etc.)
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

    %· · · · · · · · · · · · · · · · · · · · · · · · · · · ·
    % custom mouseover, mouseclick & scrollwheel functions
    %· · · · · · · · · · · · · · · · · · · · · · · · · · · ·

    % mouseover effects for Pareto front plot
    function current_Pareto = Paretos_mouseover(good_points, bad_points, varargin)
        % initially, select no pareto (same function as [dummy])
        current_Pareto = [];
        % double check if these axes are indeed the current axes
        if get(MainFig(1), 'currentaxes') ~= Pareto_axes, return; end
        % first, reset all points to normal
        set(good_points, ...
            'marker'    , '.', ...
            'color'     , [.2, .9, .0], ...
            'markersize', 15);
        set(bad_points, ...
            'marker'    , 'x', ...
            'color'     , 'r', ...
            'markersize', 5);
        % get current position
        % NOTE: HITTEST() IS AN UNDOCUMENTED FUNCTION AND MIGHT
        % DISAPPEAR WITHOUT WARNING IN A NEXT MATLAB VERSION
        close_good = (hittest == good_points);
        close_bad  = (hittest == bad_points);
        % reset the marker if the mouse hovers over it
        if any(close_good)
            % return the currently selected point
            % (make sure sure it's only 1)
            current_Pareto = good_points(find(close_good, 1));
            % set its marker
            set(current_Pareto, ...
                'marker'    , '*', ...
                'color'     , [1, .6, .2],...
                'markersize', 15);
        end
        if any(close_bad)
            % return the currently selected point
            % (make sure sure it's only 1)
            current_Pareto = bad_points(find(close_bad, 1));
            % set its marker
            set(current_Pareto, ...
                'marker'    , 'x', ...
                'color'     , [1, 0, 0],...
                'markersize', 15);
        end
        % update the info pane with the information from the point
        % the mouse currently hovers over

        %??? TODO

    end % Paretos_mouseover

    % plot new trajectory upon clicking a Pareto-solution
    function plot_new_trajectory(good_points, bad_points, varargin)
        % double check if these axes are indeed the current axes
        if get(MainFig(1), 'currentaxes') ~= Pareto_axes, return, end
        % first, get the currently "selected" pareto-point
        current_Pareto = Paretos_mouseover(good_points, bad_points);
        % no point was actually selected
        if isempty(current_Pareto), return , end
        % Emphasize it in the Pareto plot
        set(current_Pareto, ...
            'color'     , 'k', ...
            'Marker'    , '.',...
            'MarkerSize', 20);
        % we first have to find the proper index to the array
        fv1 = get(current_Pareto, 'xdata'); fv2 = get(current_Pareto, 'ydata');
        current_Pareto = ...
            results.Paretos.function_values(:, 1) == -fv1 & ...
            results.Paretos.function_values(:, 2) == fv2;
        % now plot the new trajectory associated with this new point
        results.best.solution       = results.Paretos.solutions{current_Pareto};
        results.best.function_value = results.Paretos.function_values(current_Pareto,:);
        % include new post-processor result (if applicable)
        if settings.postprocessing.check
           results.best.postprocessor_data = callback('post_processor', ...
               results.objective_function, results.Paretos.solutions{current_Pareto});
        end
        % include additional data from third objective (if applicable)
        if results.Paretos.types.other.use
            results.best.third_objective_data = ...
                results.Paretos.third_objective_data(current_Pareto,:);
        end
        % call plotting sequence
        high_first_single(true);
        % also switch to the trajectory tab
        callbacks('show_output_tab', trajectory_tab);
    end % plot_new_trajectory

    % zoom in with the mouse wheel
    % (it should be standard in MATLAB, really!)
    function dummy = scroll_zoom(axs, varargin)

        % the function should have at least one output argument
        dummy = [];

        % double check if these axes are indeed the current axes
        if get(MainFig(2), 'currentaxes') ~= axs
            return; end

        % get the amount of scolls
        scrolls = varargin{2}.VerticalScrollCount;

        % get the axes' x- and y-limits
        xlim = get(axs, 'xlim');  ylim = get(axs, 'ylim');

        % get the current camera position, and save the [z]-value
        cam_pos_Z = get(axs, 'cameraposition');  cam_pos_Z = cam_pos_Z(3);

        % get the current point
        old_position = get(axs, 'CurrentPoint'); old_position(1,3) = cam_pos_Z;

        % calculate zoom factor
        zoomfactor = max(0.01, 1 - scrolls/20);

        % adjust camera position
        set(axs,...
            'cameratarget', [old_position(1, 1:2), 0],...
            'cameraposition', old_position(1, 1:3));

        % adjust the camera view angle (equal to zooming in)
        camzoom(zoomfactor);

        % zooming with the camera has the "nasty" side-effect of
        % NOT adjusting the axes limits. We have to correct for this:
        x_lim1 = (old_position(1,1) - min(xlim))/zoomfactor;
        x_lim2 = (max(xlim) - old_position(1,1))/zoomfactor;
        xlim   = [old_position(1,1) - x_lim1, old_position(1,1) + x_lim2];
        y_lim1 = (old_position(1,2) - min(ylim))/zoomfactor;
        y_lim2 = (max(ylim) - old_position(1,2))/zoomfactor;
        ylim   = [old_position(1,2) - y_lim1, old_position(1,2) + y_lim2];
        set(axs, 'xlim', xlim), set(axs, 'ylim', ylim)

        % set new camera position
        new_position = get(axs, 'CurrentPoint');
        old_camera_target =  get(axs, 'CameraTarget');
        old_camera_target(3) = cam_pos_Z;
        new_camera_position = old_camera_target - ...
            (new_position(1,1:3) - old_camera_target(1,1:3));

        % adjust camera target and position
        set(axs, 'cameraposition', new_camera_position(1, 1:3),...
            'cameratarget', [new_camera_position(1, 1:2), 0]);

        % we also have to re-set the axes to stretch-to-fill mode
        set(axs, 'cameraviewanglemode', 'auto',...
            'camerapositionmode', 'auto',...
            'cameratargetmode', 'auto');

    end % scroll_zoom

    % pan upon mouse click
    function dummy = pan_click(axs, xlim_original, ylim_original, varargin)

        % the function should have at least one output argument
        dummy = [];

        % double check if these axes are indeed the current axes
        if get(MainWin, 'currentaxes') ~= axs
            return; end

        % perform appropriate action
        switch lower(get(MainWin, 'selectiontype'))

            % start panning on left click
            case 'normal'
                status = 'down';
                switch axs
                    case trajectory_axes
                        previous_point_trajectory = get(axs, 'CurrentPoint');
                    case Pareto_axes
                        previous_point_Paretos = get(axs, 'CurrentPoint');
                end

            % reset view on double click
            case 'open' % double click (left or right)
                set(axs, ...
                    'Xlim', xlim_original,...
                    'Ylim', ylim_original);
                axis equal

        end

    end % pan click

    % release mouse button
    function dummy = pan_release(axs, varargin)

        % assign empty dummy
        dummy = [];

        % double check if these axes are indeed the current axes
        if get(MainWin, 'currentaxes') ~= axs
            return; end

        % just reset status
        status = '';

    end % pan release

    % move the mouse (with button clicked)
    function dummy = pan_motion(axs, varargin)

        % assign empty dummy
        dummy = [];

        % double check if these axes are indeed the current axes
        if get(MainWin, 'currentaxes') ~= axs
            return; end

        % return if mouse hasn't been clicked
        if isempty(status)
            return; end

        % get current location (in pixels)
        current_point = get(axs, 'CurrentPoint');

        % get current XY-limits
        xlim = get(axs, 'xlim');  ylim = get(axs, 'ylim');

        % save new position
        switch axs
            case trajectory_axes
                % find change in position
                delta_points = current_point - previous_point_trajectory;
                % adjust limits
                new_xlim = xlim - delta_points(1);
                new_ylim = ylim - delta_points(3);
                % set new limits
                set(trajectory_axes, 'Xlim', new_xlim,...
                                    'Ylim', new_ylim);
                % save new position
                previous_point_trajectory = get(axs, 'CurrentPoint');

            case Pareto_axes
                % find change in position
                delta_points = current_point - previous_point_Paretos;
                % adjust limits
                new_xlim = xlim - delta_points(1);
                new_ylim = ylim - delta_points(3);
                % set new limits
                set(Pareto_axes, 'Xlim', new_xlim,...
                                 'Ylim', new_ylim);
                % save new position
                previous_point_Paretos = get(axs, 'CurrentPoint');
        end
    end % pan motion

end % generate output
