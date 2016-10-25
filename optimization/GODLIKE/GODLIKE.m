function varargout = GODLIKE(funfcn, lb, ub, varargin)
% GODLIKE           Global optimizer that combines the power
%                   of a Genetic Algorithm, Diffential Evolution,
%                   Particle Swarm Optimization and Adaptive
%                   Simulated Annealing number_of_algorithms.
%
% Usage:
%
% (Single-objective optimization)
%================================
%   sol = GODLIKE(obj_fun, lb, ub)
%   sol = GODLIKE(..., ub, confcn)
%   sol = GODLIKE(..., confcn, options)
%   sol = GODLIKE(..., confcn, 'option', value, ...)
%
%   [sol, fval] = GODLIKE(...)
%   [sol, fval, exitflag] = GODLIKE(...)
%   [sol, fval, exitflag, output] = GODLIKE(...)
%
%
% (Multi-objective optimization)
% ==============================
%   sol = GODLIKE(obj_fun12..., lb, ub)
%   sol = GODLIKE({obj_fun1, obj_fun2,...}, lb, ub)
%   sol = GODLIKE(..., ub, confcn, options)
%   sol = GODLIKE(..., ub, confcn, 'option', value, ...)
%
%   [sol, fval] = GODLIKE(...)
%   [..., fval, Pareto_front] = GODLIKE(...)
%   [..., Pareto_front, Pareto_Fvals] = GODLIKE(...)
%   [..., Pareto_Fvals, exitflag] = GODLIKE(...)
%   [..., exitflag, output] = GODLIKE(...)
%
%
% INPUT ARGUMENTS:
% ================
%
%   obj_fun     The objective function of which the global minimum
%               will be determined (function_handle). For multi-
%               objective optimization, several objective functions
%               may be provided as a cell array of function handles,
%               or alternatively, in a single function that returns
%               the different function values along the second
%               dimension.
%               Objective functions must accept either a [popsize x
%               dimensions] matrix argument, or a [1 x dimensions]
%               vector argument, and return a [popsize x number of
%               objectives] matrix or [1 x number of objective]
%               vector of associated function values (number of
%               objectives may be 1). With the first format, the
%               function is evaluated vectorized, in  the second
%               case CELLFUN() is used, which is a bit slower in
%               general.
%
%   lb, ub      The lower and upper bounds of the problem's search
%               space, for each dimension. May be scalar in case all
%               bounds in all dimensions are equal. Note that at
%               least ONE of these must have a size of [1 x
%               dimensions], since the problem's dimensionality is
%               derived from it.
%
%     confcn   Constraint functions.
%
%objective function of which the global minimum
%               will be determined (function_handle). For multi-
%               objective optimization, several objective functions
%               may be provided as a cell array of function handles,
%               or alternatively, in a single function that returns
%               the different function values along the second
%               dimension.
%               Objective functions must accept either a [popsize x
%               dimensions] matrix argument, or a [1 x dimensions]
%               vector argument, and return a [popsize x number of
%               objectives] matrix or [1 x number of objective]
%               vector of associated function values (number of
%               objectives may be 1). With the first format, the
%               function is evaluated vectorized, in  the second
%               case CELLFUN() is used, which is a bit slower in
%               general.
%
%   which_ones  The number_of_algorithms to be used in the optimizations. May
%               be a single string, e.g., 'DE', in which case the
%               optimization is equal to just running a single
%               Differential Evolution optimization. May also be a
%               cell array of strings, e.g., {'DE'; 'GA'; 'ASA'},
%               which uses all the indicated number_of_algorithms. When
%               omitted or left empty, defaults to {'DE';'GA';'PSO';
%               'ASA'} (all number_of_algorithms once).
%
%   options/    Sets the options to be used by GODLIKE. Options may
%   'option',   be a structure created by set_options, or given as
%      value    individual ['option', value] pairs. See set_options
%               for a list of all the available options and their
%               defaults.
%
% OUTPUT ARGUMENTS:
% =================
%
%   sol         The solution that minizes the problem globally,
%               of size [1 x dimensions]. For multi-objective
%               optimization, this indicates the point with the
%               smallest distance to the (shifted) origin.
%
%   fval        The globally minimal function value
%
%   exitflag    Additional information to facilitate fully automated
%               optimization. Negative is `bad', positive `good'. A
%               value of '0' indicates GODLIKE did not perform any
%               operations and exited prematurely. A value of '1'
%               indicates normal exit conditions. A value of '2' meand
%               that one of the provided output functions issued a
%               stop-command. A value of '-1' indicates a premature
%               exit due to exceeding the preset maximum number of
%               function evaluations. A value of '-2' indicates that
%               the amount of maximum GODLIKE iterations has been
%               exceeded, and a value of '-3' indicates no optimum
%               has been found (only for single-objective
%               optimization).
%
%   output      structure, containing much additional information
%               about the optimization as a whole; see the manual
%               for a more detailed description.
%
%   (For multi-objective optimization only)
%
%   Pareto_front, Pareto_Fvals
%               The full set of non-dominated solutions, and their
%               associated function values.
%
%   See also pop_single, pop_multi, set_options.


%      Author : Rody P.S. Oldenhuis
% Affiliation : Delft University of Technology
%               Faculty of Aerospace Engineering
%               Dep. of Astrodynamics & Satellite Systems
%     Contact : oldenhuis@gmail.com
%   Licensing/
%    (C) info : Frankly I don't care what you do with it,
%               as long as I get some credit when you copy
%               large portions of the code ^_^

% last edited 03/10/2009

    %% INITIALIZE

    % basic check on input
    error(nargchk(3, inf, nargin));

    % more elaborate check on input (nested function)
    check_input;

    % resize and reshape boundaries and dimensions (nested function)
    [lb, ub, sze, popsize, dimensions, confcn, constrained, which_ones, options] = ...
        reformat_input(lb, ub, varargin{:});

    % test input objective function(s) to determine the problem's dimensions,
    % number of objectives and proper input format (nested function)
    [options, single, multi, test_evaluations] = test_funfcn(options);

    % initialize more variables
    number_of_algorithms = numel(which_ones);            % number of number_of_algorithms to use
       number_of_streams = options.NumStreams;           % number of simultaneous streams
              generation = 1;                            % this is the first generation
     num_funevaluations  = 0;                            % number of function evaluations
     [converged, output] = check_convergence;            % initial output structure
          outputFcnbreak = false;                        % exit condition for output functions
                     pop = cell(number_of_streams, number_of_algorithms);
                                                         % cell array of population-objects

    % Initially, [output] is a large structure used to move data to and from all the
    % subfunctions. Later, it is turned into the output argument [output] by removing some
    % obsolete entries from the structure.

    % if an output function's been given, evaluate them
    state = 'init'; % initialization state
    if ~isempty(options.outputFcn)
        cellfun(@(x) x([],[],state), options.outputFcn, 'uniformoutput', false);
    end

    % do an even more elaborate check (the behavior of this
    % nested function is changed by passing the number of
    % requested output arguments)
    check_input(nargout);

    %% GODLIKE loop

    % GODLIKE loop
    while ~converged

        % randomize population sizes for each stream
        % (minimum is 5 individuals)
        frac_popsize = break_value(popsize, number_of_streams, 5);

        % shuffle populations (NOT in the first iteration)
        if (generation > 1), pop = interchange_populations(pop); end

        % loop through each stream
        for i = 1:number_of_streams

            % randomize number of iterations per algorithm
            % ([options.ItersUb] is the maximum TOTAL amount
            % of iterations that will be spent in all of the algorithms combined.)
            frac_iterations = ...
                break_value(options.ItersUb, number_of_algorithms, options.ItersLb);

            % loop through each algorithm
            for j = 1:number_of_algorithms

                % initialize this population (if required)
                pop = initialize_population(pop);

                % Perform single iterations
                counter = 0; % used for single-objective optimization
                for k = 1:frac_iterations(j)

                    % do single iteration on this population
                    pop{i,j}.iterate;

                    % evaluate the output functions
                    if ~isempty(options.outputFcn)
                        % most intensive part, here in the inner loop
                        state = 'interrupt';
                        % collect information
                        [x, optimValues] = get_outputFcn_values(i,j);
                        % evaluate the output functions
                        stop = cellfun(@(y)y(x, optimValues, state), ...
                            options.outputFcn, 'uniformoutput', false);
                        stop = any([stop{:}]);
                        % GODLIKE might need to stop here
                        if stop, outputFcnbreak = true; break; end
                    end

                    % calculate total number of function evaluations
                    initialized_populations = ~cellfun(@isempty, pop);
                    funevaluations = cellfun(@(x)x.funevals, pop(initialized_populations));
                    num_funevaluations = test_evaluations + sum(funevaluations);

                    % check for convergence of this iteration
                    if multi
                        % all are non-dominated, first display progress, then exit the loop
                        [alg_converged, output] = check_convergence(false,output,[]);
                        if alg_converged
                            if ~isempty(options.display), display_progress; end
                            break
                        end
                    elseif single
                        % check algorithm convergence
                        [alg_converged, output, counter] = check_convergence(false,output,counter);
                        % if converged, first display progress, empty the other
                        % populations in this stream, and exit the loop
                        if alg_converged
                            if ~isempty(options.display), display_progress; end
                            pop(i, (j+1):end) = {[]};
                            break
                        end
                    end % if

                    % check function evaluations, and exit if it
                    % surpasses the preset maximum
                    if (num_funevaluations >= options.MaxFunEvals)
                        % also display last iteration
                        if ~isempty(options.display), display_progress; end,
                        converged = true; break
                    end

                    % display progress at every iteration
                    if ~isempty(options.display), display_progress; end

                end % algorithm inner loop

                % if one of the output functions returned a stop request, break
                if outputFcnbreak, break, end
                % if we have convergence inside the algorithm inner loop, break
                if alg_converged, break; end
                % if the max. # of function evaluations has been surpassed, break
                if converged, break; end

                % evaluate the output functions
                if ~isempty(options.outputFcn)
                    % end of an algorithm loop
                    state = 'iter';
                    % collect the information
                    [x, optimValues] = get_outputFcn_values(i,j);
                    % call the output functions
                    cellfun(@(y)y(x,optimValues,state), options.outputFcn, 'uniformoutput', false);
                end

            end % algorithm outer loop

            % if one of the output functions returned a stop request, break
            if outputFcnbreak, break, end
            % if the max. # of function evaluations has been surpassed, break
            if converged, break; end
            % if we have convergence inside the algorithm loop,
            % determine whether we should quit or not
            if alg_converged
                % don't continue if we only have one stream
                if number_of_streams == 1, converged = true; break; end
                % don't continue if QuiWhenAchieved is on, and the value
                % has been achieved
                if options.QuitWhenAchieved && single && ...
                   output.global_best_funval <= options.AchieveFunVal
                        converged = true; break
                end
                % DO continue if these conditions do not hold
                alg_converged = false;
            end

        end % stream loop

        % if one of the output functions returned a stop request, break
        if outputFcnbreak, converged = true; end
        % check maximum iterations
        if (generation >= options.MaxIters), converged = true; end

        % check for GODLIKE convergence (and update output structure)
        [converged, output] = check_convergence(converged, output);

        % increase generation
        generation = generation + 1;

        % evaluate the output functions
        if ~outputFcnbreak && ~isempty(options.outputFcn)
            % end of a GODLIKE iteration
            state = 'iter';
            % collect the information
            [x, optimValues] = get_outputFcn_values([],[]);
            % call the output functions
            cellfun(@(y)y(x, optimValues, state), options.outputFcn, 'uniformoutput', false);
        end

    end % GODLIKE loop

    % display final results
    % (*NOT* if the output function requested to stop)
    if ~outputFcnbreak && ~isempty(options.display), display_progress; end

    %% OUTPUT VALUES

    % multi-objective optimization
    if multi
        varargout{1} = reshape(output.most_efficient_point, sze);
        varargout{2} = output.most_efficient_fitnesses;
        varargout{3} = reshape(output.pareto_front_individuals.', ...
            [sze, size(output.pareto_front_individuals, 1)]);
        varargout{4} = output.pareto_front_fitnesses;
        varargout{5} = output.exitflag;
        % remove some fields from output structure
        output = rmfield(output, {'pareto_front_individuals','pareto_front_fitnesses',...
            'exitflag','most_efficient_point','most_efficient_fitnesses'});
        % and output what's left
        varargout{6} = output;
        % in case the output function requested to stop
        if outputFcnbreak
            output.exitflag = 2;
            output.message  = 'GODLIKE was terminated by one of the output functions.';
            varargout{5} = output.exitflag ;
            varargout{6} = output;
        end

    % single-objective optimization
    elseif single

        % if all went normal
        if isfield(output, 'global_best_individual')
            varargout{1} = reshape(output.global_best_individual, sze);
            if constrained
                varargout{2} = output.global_best_funval_unconstrained;
            else
                varargout{2} = output.global_best_funval;
            end
            varargout{3} = output.exitflag;
            % remove some fields from output structure
            outpt.algorithms = output.algorithms;   outpt.funcCount = output.funcCount;
            outpt.message    = output.message;      outpt.algorithm_info = output.algorithm_info;
            outpt.iterations = output.iterations;
            % and output
            varargout{4} = outpt;
            % in case the output function requested to stop
            if outputFcnbreak
                outpt.exitflag = 2;
                outpt.message  = 'GODLIKE was terminated by one of the output functions.';
                varargout{3} = outpt.exitflag;
                varargout{4} = outpt;
            end

        % but, no optimum might have been found
        else
            varargout{1} = NaN(sze);
            varargout{2} = NaN;
            varargout{3} = -3;
            % remove some fields from output structure
            output = rmfield(output, {'global_best_funval', 'exitflag','descent_counter',...
                'best_individuals','best_funcvalues','previous_global_best_funval',...
                'previous_best_funcvalues'});
            % adjust message
            output.message = sprintf('%s\n\n All function values encountered were INF or NaN.\n',...
                output.message);
            % output
            varargout{4} = output;
            % in case the output function requested to stop
            if outputFcnbreak
                output.message  = 'GODLIKE was terminated by one of the output functions.';
                output.exitflag = 2;
                varargout{3} = output.exitflag;
                varargout{4} = output;
            end
        end
    end

    % last call to output function
    if ~isempty(options.outputFcn)
        cellfun(@(y)y([],[], 'done'), options.outputFcn, 'uniformoutput', false);
    end

    %% NESTED FUNCTIONS

    % =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =
    % initialization shizzle
    % =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    % elaborate error trapping
    function check_input(varargin)
        if (nargin == 0)
            if isempty(funfcn)
                error('GODLIKE:function_not_defined',...
                    'GODLIKE requires at least one objective function.');
            end
            if isempty(lb)||isempty(ub)
                error('GODLIKE:lbubpopsize_not_defined',...
                    'GODLIKE requires arguments [lb] and [ub].');
            end
            if ~isnumeric(lb)||~isnumeric(ub)
                error('GODLIKE:lbubpopsize_not_numeric',...
                    'Arguments [lb], and [ub] must be numeric.');
            end
            if any(~isfinite(lb)) || any(~isfinite(ub)) || ...
                    any(  ~isreal(lb)) || any(~isreal(ub))
                error('GODLIKE:lbub_not_finite',...
                    'Values for [lb] and [ub] must be real and finite.');
            end
            if ~isa(funfcn, 'function_handle')
                % might be cell array
                if iscell(funfcn)
                    for ii = 1:numel(funfcn)
                        if ~isa(funfcn{ii}, 'function_handle')
                            error('GODLIKE:funfcn_mustbe_function_handle',...
                                'All objective functions must be function handles.');
                        end
                    end
                % otherwise, must be function handle
                else
                    error('GODLIKE:funfcn_mustbe_function_handle',...
                        'Objective function must be given as a function handle.');
                end
            end
            if (nargin == 7) && ~isstruct(varargin{2})
                error('GODLIKE:options_mustbe_structure',...
                    'Argument [options] must be a structure.')
            end
            if any(any(lb > ub))
                error('GODLIKE:lb_larger_than_ub',...
                    'All entries in [lb] must be smaller than the corresponding entries in [ub].')
            end
        else
            if (options.ItersLb > options.ItersUb)
                warning('GODLIKE:ItersLb_exceeds_ItersUb',...
                    ['Value of options.ItersLb is larger than value of\n',...
                    'options.ItersUb. Values will simply be swapped.']);
                u_b = options.ItersUb;
                options.ItersUb = options.ItersLb;
                options.ItersLb = u_b;
            end
            if (options.ItersLb > options.ItersUb)
                warning('GODLIKE:MaxIters_exceeds_MinIters',...
                    ['Value of options.MinIters is larger than value of\n',...
                    'options.MaxIters. Values will simply be swapped.']);
                u_b = options.MaxIters;
                options.MaxIters = options.MinIters;
                options.MinIters = u_b;
            end
            if single
                % single objective optimization has a maximum of 4 output arguments
                error(nargoutchk(0, 4, varargin{1}))
            elseif multi
                % multi-objective optimization has a maximum of 6 output arguments
                error(nargoutchk(0, 6, varargin{1}))
            end
            if ~isempty(options.outputFcn) && ...
               ~all( cellfun(@(x) isa(x, 'function_handle'), options.outputFcn))
                error('GODLIKE:outputFcn_shouldbe_function_handle',...
                    'All output functions should be function handles.')
            end
        end % if
    end % nested function

    % reshape, resize and redefine input to predictable formats
    function [lb, ub, sze, popsize, dimensions, confcn, constrained, which_ones, options] = ...
             reformat_input(lb, ub, varargin)

        % First set /get options
        if nargin <= 3, options = set_options; end                  % defaults
        if nargin == 4, options = varargin{2}; end                  % structure provided
        if nargin > 4 , options = set_options(varargin{2:end}); end % individually provided

        % cast output functions to cell
        if isfield(options, 'outputFcn') && isa(options.outputFcn, 'function_handle')
            options.outputFcn = {options.outputFcn}; end

        % constraint functions
        if nargin == 2 || isempty(varargin{1}) % default - no constraint function
            confcn = {[]};
            % constraint might also be calculated inside the objective
            % function(s)
            if options.ConstraintsInObjectiveFunction > 0, constrained = true;
            else constrained = false;
            end
        else
            confcn = varargin{1};
            constrained = true;
            % cast to cell if only one is selected
            if isa(confcn, 'function_handle'), confcn = {confcn}; end
            % possible erroneous input
            if ~iscell(confcn)
                error('GODLIKE:confcn_mustbe_cell_or_funhandle',...
                    ['Constraint functions must be given as a fuction_handle, or as a cell ',...
                    'array of function_handles.']);
            end
            % FURTHER CHECKING WILL BE DONE IN TEST_FUNFCN()
        end

        % extract which algorithms to use
        which_ones = options.algorithms;

        % save the original size of [lb] or [ub]
        max_elements = max(numel(lb),numel(ub));
        if (max_elements == numel(lb)), sze = size(lb); else  sze = size(ub); end

        % force [lb] and [ub] to be row vectors
        lb = lb(:).';  ub = ub(:).';

        % replicate one or the other when their sizes are not equal
        if ~all(size(lb) == size(ub))
            if     isscalar(lb)
                lb = repmat(lb, size(ub));
            elseif isscalar(ub)
                ub = repmat(ub, size(lb));
            else
                error('GODLIKE:lbub_sizes_incorrect',...
                     ['If the size of either [lb] or [ub] is equal to the problem''s dimenions\n',...
                     'the size of the other must be 1x1.'])
            end
        end

        % define [dimensions]
        dimensions = numel(lb);

        % total population size
        % (defaults to 25*number of dimensions)
        if isempty(options.popsize)
            popsize = min(25*dimensions, 1500);
        else
            popsize = options.popsize;
        end

        % check minimum popsize
        minpop = 5*numel(options.algorithms);
        if (minpop > popsize)
            warning('GOLIKE:popsize_too_small',...
                ['Each algorithm requires a population size of at least 5.\n',...
                'Given value for [popsize] makes this impossible. Increasing\n',...
                'argument [popsize] to ', num2str(minpop), '...']);
            popsize = minpop;
        end

        % Assign back for consistency
        options.GODLIKE.popsize = popsize;

    end % reformat input

    % test the function, and determine the amount of objectives. Here
    % it is decided whether the optimization is single-objective or
    % multi-objective. Also test the given constraint function
    function [options, single, multi, fevals] = test_funfcn(options)

        % initialize
        fevals = 0; options.num_objectives = 1;

        % split multi/single objective
        if iscell(funfcn) && (numel(funfcn) > 1)
            % no. of objectives is simply the amount of provided objective functions
            options.num_objectives = numel(funfcn);
            % single is definitely false
            single = false;
        elseif iscell(funfcn) && (numel(funfcn) == 1)
            % single it true but might still change to false
            single = true;
        else
            % cast [funfcn] to cell
            funfcn = {funfcn};
            % single is true but might still change to false
            single = true;
        end

        % � � � � � � � � � � � � � � � � � � � � � � � � � � � � � � � �
        % try to evaluate the objective and constraint functions, with the
        % original [lb]. If any evaluation fails, throw an error
        % � � � � � � � � � � � � � � � � � � � � � � � � � � � � � � � �

        % reshape to original size
        lb_original = reshape(lb, sze);

        % loop through (all) objective function(s)
        for ii = 1:numel(funfcn)

            % try to evaluate the function
            try
                % simply evaluate the function with the lower bound
                if options.ConstraintsInObjectiveFunction == 0
                    % no constraints
                    sol = feval(funfcn{ii}, (lb_original));
                % constraints might be given in the objective functions
                else
                    arg_out = cell(1, options.ConstraintsInObjectiveFunction);
                    [arg_out{:}] = feval(funfcn{ii}, lb_original);
                    sol = arg_out{1};
                    con = arg_out{options.ConstraintsInObjectiveFunction};
                    % con MUST be a vector
                    if ~isvector(con)
                        error('GODLIKE:confun_must_return_vector',...
                            ['All constraint functions must return a [Nx1] or [1xN] vector ',...
                            'of constraint violations.\nSee the documentation for more details.']);
                    end
                end
                % keep track of the number of function evaluations
                fevals = fevals + 1;

                % see whether single must be changed to multi
                if single && (numel(sol) > 1)
                    single = false;
                    options.num_objectives = numel(sol);
                end

                % it might happen that more than one function is provided,
                % but that one of the functions returns more than one function
                % value. GODLIKE does not handle that case
                if (numel(sol) > 1) && (ii > 1)
                    error('GODLIKE:multimulti_not_allowed',...
                        ['GODLIKE cannot optimize multiple multi-objective problems ',...
                        'simultaneously.\nUse GODLIKE multiple times on each of your objective ',...
                        'functions separately.\n\nThis error is generated because the first of ',...
                        'your objective functions returned\nmultiple values, while ',...
                        'you provided multiple objective functions. Only one of\nthese formats ',...
                        'can be used for multi-objective optimization, not both.'])
                end

            % if evaluating the function fails, throw error
            catch userFcn_ME
                pop_ME = MException('GODLIKE:function_doesnt_evaluate',...
                    'GODLIKE cannot continue: failure during function evaluation.');
                userFcn_ME = addCause(userFcn_ME, pop_ME);
                rethrow(userFcn_ME);
            end % try/catch

        end % for

        % see if the optimization is multi-objective
        multi = ~single;

        % loop through (all) constraint function(s)
        if constrained && (options.ConstraintsInObjectiveFunction == 0)
            for ii = 1:numel(confcn)
                % some might be empty
                if isempty(confcn{ii}), continue, end
                % try to evaluate the function
                try
                    % simply evaluate the function with the lower bound
                    con = feval(confcn{ii}, lb_original);
                    % keep track of the number of function evaluations
                    fevals = fevals + 1;
                    % con MUST be a vector
                    if ~isvector(con)
                        error('GODLIKE:confun_must_return_vector',...
                            ['All constraint functions must return a [Nx1] or [1xN] vector ',...
                            'of constraint violations.\nSee the documentation for more details.']);
                    end

                % if evaluating the function fails, throw error
                catch userFcn_ME
                    pop_ME = MException('GODLIKE:constraint_function_doesnt_evaluate',...
                        'GODLIKE cannot continue: failure during evaluation of one of the constraint functions.');
                    userFcn_ME = addCause(userFcn_ME, pop_ME);
                    rethrow(userFcn_ME);
                end % try/catch
            end % for
        end % if constrained

    end % test function

    % =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =
    % functions used in the main loop
    % =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    % break up some [value] into a vector of random integers
    % of length [number_of_algorithms], that sums up to [value]
    function frac_value = break_value(value, number_of_values, Lb)
        % NOTE: The case of these variables [Lb] and [Ub] is important.
        % The GODLIKE arguments [lb] or [ub] may get overwritten!

        % only one algorithm - just return value
        if number_of_values == 1, frac_value = value; return; end

        % initially, the upper bound is the value minus
        % (number_of_values-1) times the lower bound
        Ub = value - (number_of_values-1)*Lb;

        % create array of length [number_of_algorithms] that
        % sums to [value]
        frac_value = zeros(number_of_values, 1);
        for ii = 1:number_of_values-1 % note the minus one
            % random value (make sure it's not zero)
            rnd = 0; while (rnd == 0), rnd = round(rand*(Ub-Lb) + Lb); end
            frac_value(ii) = rnd;
            % adjust max. value for next iteration
            Ub = round((value - sum(frac_value))/(number_of_values-ii));
        end % for

        % last entry is the difference of the sum of all values and the original value
        frac_value(end) = value - sum(frac_value);

        % sort at random
        [dummy, inds] = sort(rand(size(frac_value,1),1));%#ok
        frac_value = frac_value(inds);

    end % break value

    % initialize populations
    function pop = initialize_population(pop)
        % rename for clarity
        stream = i;  algorithm = j;

        % initialize first populations (if this is the first call)
        if (generation == 1) && (stream == 1) && (algorithm == 1) && isempty(pop{1,1})
            % set algorithm
            type = upper(which_ones{1});
            options.algorithm = type;
            % initialize all first populations
            for ii = 1:number_of_streams
%                 [dummy, pop{ii, 1}] = ...
%                     evalc([type, '(funfcn,confcn,frac_popsize(ii),lb,ub,sze,options)']);
pop{ii, 1} = eval([type, '(funfcn,confcn,frac_popsize(ii),lb,ub,sze,options)']);
            end
            return % we're done

        % if this is not the first call, pass results from the
        % previous population into the next one
        elseif (algorithm > 1)
            % set algorithm
            type = upper(which_ones{algorithm});
            options.algorithm = type;
            % previous algorithm
            previous_algorithm = algorithm - 1;
            % initialize the next algorithm
            prev_pop = pop{stream,previous_algorithm};%#ok
            [dummy, pop{stream, algorithm}] = ...
                evalc([type, '(prev_pop.pop_data,prev_pop,options)']);%#ok

        % exceptional case: at a next generation, the function count needs
        % to be updated, but NOT the population (as this is already done in
        % INTERCHANGE_POPULATIONS())
        elseif (generation > 1) && (algorithm == 1)
            % find the last non-empty population
            last_nonempty_pop = sum(~cellfun(@isempty, pop(stream, :)));
            % THIS is the last algorithm
            previous_algorithm = last_nonempty_pop;
            % get the previous function count
            prev_pop_funcCount = pop{stream,previous_algorithm}.funevals;
            % set the next algorithm's function count
            pop{stream,1}.funevals = prev_pop_funcCount;
        end
    end % initialize population

    % shuffle and (re)initialize the population objects
    function pop = interchange_populations(pop)

        % rename for clarity
        stream = i;

        % if there's only one stream, copy the last non-empty algorithms
        % to the first one, and return (we don't need to INTERCHANGE in
        % this case)
        if (number_of_streams == 1)
            % find the last non-empty population
            last_nonempty_pop = sum(~cellfun(@isempty, pop(stream, :)));
            % get the type of the first algorithm in the stream
            type = upper(which_ones{1});
            % THIS is the last algorithm
            previous_algorithm = last_nonempty_pop;
            % get its population
            prev_pop = pop{1,previous_algorithm};%#ok
            % initialize the first algorithm
            [dummy, pop{1,1}] = ...
                evalc([type, '(prev_pop.pop_data,prev_pop,prev_pop.options)']);%#ok
            return % we're done
        end

        % initialize
        parent_pops    = zeros(popsize, pop{1}.dimensions);
        parent_fits    = zeros(popsize, options.num_objectives);
        offspring_pops = parent_pops;
        offspring_fits = parent_fits;
        if multi
            front_numbers      = zeros(popsize, 1);
            crowding_distances = [front_numbers;front_numbers];
        end
        if constrained
            parent_constrviolation     = zeros(popsize, numel(confcn));
            parent_unpenalized_fits    = zeros(popsize, options.num_objectives);
            offspring_constrviolation  = parent_constrviolation;
            offspring_unpenalized_fits = parent_unpenalized_fits;
        end
        lfe1 = 0;   lfe2 = 0;    % Last Filled Entry (lfe)

        % extract all current populations, their function values,
        % and other relevant information
        for ii = 1:number_of_streams

            % find the last non-empty population in this stream
            jj = sum(~cellfun(@isempty, pop(ii, :)));

            % rename stuff for clarity
            popinfo = pop{ii,jj}.pop_data;       popsz = pop{ii,jj}.size;

            % both for single and multi-objective
            parent_pops(lfe1+1:lfe1+popsz, :)  = popinfo.parent_population;
            parent_fits(lfe1+1:lfe1+popsz, :)  = popinfo.function_values_parent;
            offspring_pops(lfe1+1:lfe1+popsz,:)= popinfo.offspring_population;
            offspring_fits(lfe1+1:lfe1+popsz,:)= popinfo.function_values_offspring;

            % stuff specific for multi-objective optimization
            if multi
                front_numbers(lfe1+1:lfe1+popsz, :)        = popinfo.front_number;
                crowding_distances(lfe2+1:lfe2+2*popsz, :) = popinfo.crowding_distance;
            end

            % stuff specific for constrained optimization
            if constrained
                parent_constrviolation(lfe1+1:lfe1+popsz, :)     = popinfo.constraint_violations_parent;
                parent_unpenalized_fits(lfe1+1:lfe1+popsz, :)    = popinfo.unpenalized_function_values_parent;
                offspring_constrviolation(lfe1+1:lfe1+popsz, :)  = popinfo.constraint_violations_offspring;
                offspring_unpenalized_fits(lfe1+1:lfe1+popsz, :) = popinfo.unpenalized_function_values_offspring;
            end

            % update indices
            lfe1 = lfe1 + popsz;  lfe2 = lfe2 + 2*popsz;

        end % for

        % re-initialize populations accordingly
        for ii = 1:number_of_streams

            % rename for clarity
            fp = frac_popsize(ii);

            % split everything up according to current [frac_popsize]
            new_popinfo.parent_population         = parent_pops(1:fp, :);
            new_popinfo.function_values_parent    = parent_fits(1:fp, :);
            new_popinfo.offspring_population      = offspring_pops(1:fp, :);
            new_popinfo.function_values_offspring = offspring_fits(1:fp, :);
            if multi
                new_popinfo.front_number      = front_numbers(1:fp, :);
                new_popinfo.crowding_distance = crowding_distances(1:2*fp, :);
            end % if

            % stuff specific for constrained optimization
            if constrained
                new_popinfo.constraint_violations_parent = ...
                    parent_constrviolation(1:fp, :);
                new_popinfo.unpenalized_function_values_parent = ...
                    parent_unpenalized_fits(1:fp, :);
                new_popinfo.constraint_violations_offspring = ...
                    offspring_constrviolation(1:fp, :);
                new_popinfo.unpenalized_function_values_offspring = ...
                    offspring_unpenalized_fits(1:fp, :);
            end % if constrained

            % change options - options for ASA are always different
            opts = pop{ii,1}.options;
            % apply re-heating
            opts.ASA.T0 = opts.ASA.T0 / options.ASA.ReHeating;
            % re-initialize
            type = which_ones{1};
            [dummy, pop{ii,1}] = ...
                evalc([type, '(new_popinfo, pop{ii,1}, opts)']);%#ok

            % shrink arrays (using "... = [];" for deletion is rather slow)
            parent_pops = parent_pops(fp+1:end,:);  offspring_pops = offspring_pops(fp+1:end,:);
            parent_fits = parent_fits(fp+1:end,:);  offspring_fits = offspring_fits(fp+1:end,:);
            if multi
                front_numbers      = front_numbers(fp+1:end,:);
                crowding_distances = crowding_distances(fp+1:end,:);
            end % if

        end % for

    end % nested function

    % update output values, and check for convergence
    function [converged, output, counter] = ...
             check_convergence(converged, output, varargin)

         % some number_of_algorithms might be doubly used.
         % save which ones they are
         persistent sames

         % no input - initialize
         if (nargin == 0)

             % initially, no convergence
             converged = false;

             % some number_of_algorithms might be doubly used. Find out
             % which ones, and create proper indices
             sames = ones(number_of_algorithms, 1);
             for ii = 1:number_of_algorithms
                 same        = strcmpi(which_ones, which_ones{ii});
                 sames(same) = 1:nnz(same);
             end

             % general settings
             output.algorithms = upper(which_ones); % algorithms used
             output.exitflag   = 0;                 % neutral exitflag
             output.message    = sprintf('No iterations have been performed.');
             output.funcCount  = 0;
             for ii = 1:number_of_streams
                 for jj = 1:number_of_algorithms
                     if number_of_streams > 1
                         output.algorithm_info.(...
                             ['Stream',num2str(ii)]).(upper(which_ones{jj}))(...
                             sames(jj)).funcCount  = 0;
                         output.algorithm_info.(...
                             ['Stream',num2str(ii)]).(upper(which_ones{jj}))(...
                             sames(jj)).iterations = 0;
                     else
                         output.algorithm_info.(upper(which_ones{ii}))(sames(jj)).funcCount  = 0;
                         output.algorithm_info.(upper(which_ones{ii}))(sames(jj)).iterations = 0;
                     end
                 end
             end

             % initialize [output] for single-objective optimization
             if single
                 output.descent_counter               = 0;
                 output.global_best_individual        = NaN(1,dimensions);
                 output.previous_global_best_individual = NaN(1,dimensions);
                 output.global_best_funval            = inf;
                 output.previous_global_best_funval   = inf;
                 output.best_funcvalues               = inf(number_of_streams,number_of_algorithms);
                 output.previous_best_funcvalues      = inf(number_of_streams,number_of_algorithms);
                 output.best_individuals              = NaN(number_of_streams,number_of_algorithms,dimensions);
                 output.previous_best_individuals     = NaN(number_of_streams,number_of_algorithms,dimensions);
                 % constrained problems
                 if constrained
                     output.global_best_funval_unconstrained   = inf;
                     output.global_best_funval_constrviolation = inf;
                 end
                 for ii = 1:number_of_streams
                     for jj = 1:number_of_algorithms
                         if number_of_streams > 1
                             output.algorithm_info.(...
                                 ['Stream',num2str(ii)]).(upper(which_ones{jj}))(...
                                 sames(jj)).last_population = [];
                             output.algorithm_info.(...
                                 ['Stream',num2str(ii)]).(upper(which_ones{jj}))(...
                                 sames(jj)).last_fitnesses  = [];
                         else
                             output.algorithm_info.(upper(which_ones{jj}))(...
                                 sames(jj)).last_population = [];
                             output.algorithm_info.(upper(which_ones{jj}))(...
                                 sames(jj)).last_fitnesses  = [];
                         end
                     end
                 end
             end

             % initialize [output] for multi-objective optimization
             if multi
                 output.pareto_front_individuals = [];
                 output.pareto_front_fitnesses   = [];
                 output.most_efficient_point     = [];
                 output.most_efficient_fitnesses = [];
             end

             % we're done
             return
         end

         % both per-algorithm and global check needs to be performed.
         % the mode of operation depends on the presence of a third
         % input argument. If given, only the current populations is
         % checked. If omitted, all populations are checked.
         stream = i;  algorithm = j;
         if (nargin == 3), alg_conv = true; counter  = varargin{1};
         else alg_conv = false;
         end

         % adjust iterations etc.
         output.funcCount  = num_funevaluations;
         output.iterations = generation;
         for ii = 1:number_of_streams
             for jj = 1:number_of_algorithms
                 if ~isempty(pop{ii,jj})
                     if number_of_streams > 1
                         output.algorithm_info.(['Stream',...
                             num2str(ii)]).(upper(which_ones{jj}))(sames(jj)).iterations = ...
                             pop{ii,jj}.iterations;
                         output.algorithm_info.(['Stream',...
                             num2str(ii)]).(upper(which_ones{jj}))(sames(jj)).funcCount = ...
                             pop{ii,jj}.funevals;
                     else
                         output.algorithm_info.(...
                             upper(which_ones{ii}))(sames(jj)).iterations = pop{ii,jj}.iterations;
                         output.algorithm_info.(...
                             upper(which_ones{ii}))(sames(jj)).funcCount  = pop{ii,jj}.funevals;
                     end
                 end
             end
         end

         % convergence might already have occured. Determine the reason
         % (this is both for multi- and single-objective optimizations)
         if converged
             % maximum function evaluations have been exceeded.
             if (num_funevaluations >= options.MaxFunEvals)
                 output.exitflag = -1;
                 output.message = sprintf(['Optimization terminated:\n',...
                     ' Maximum amount of function evaluations has been reached.\n',...
                     ' Increase ''MaxFunEvals'' option.']);
             end % if
             % maximum allowable iterations have been exceeded.
             if (generation >= options.MaxIters)
                 output.exitflag = -2;
                 output.message = sprintf(['Optimization terminated:\n',...
                     ' Maximum amount of iterations has been reached.\n',...
                     ' Increase ''MaxIters'' option.']);
             end % if
         end % if

         % stuff specific for single objective optimization
         if single

             % store previous global best function value
             output.previous_global_best_individual = output.global_best_individual;
             output.previous_global_best_funval     = output.global_best_funval;
             output.previous_best_funcvalues        = output.best_funcvalues;
             output.previous_best_individuals       = output.best_individuals;

             % assign global best individuals and their function
             % values per algorithm
             ind = ones(number_of_streams,number_of_algorithms);
             for ii = 1:number_of_streams
                 % find the last non-empty population in this stream
                 jj = sum(~cellfun(@isempty, pop(ii, :)));
                 % adjust the output
                 [output.best_funcvalues(ii,jj), ind(ii,jj)] = min(pop{ii,jj}.fitnesses);
                 output.best_individuals(ii,jj,:) = pop{ii,jj}.individuals(ind(ii,jj), :);
             end

             % save new global best individual and function value
             [min_func_val, index_x] = min(output.best_funcvalues,[],1);
             [min_func_val, index_y] = min(min_func_val,[],2);
             index_x = index_x(index_y);
             if (min_func_val < output.global_best_funval)
                 output.global_best_funval     = min_func_val;
                 output.global_best_individual = output.best_individuals(index_x,index_y, :);
                 % also save completely unpenalized values for constrained functions
                 if constrained
                     output.global_best_funval_unconstrained = ...
                         pop{index_x,index_y}.pop_data.unpenalized_function_values_parent(ind(index_x,index_y));
                     output.global_best_funval_constrviolation = ...
                         pop{index_x,index_y}.pop_data.constraint_violations_parent(ind(index_x,index_y), :);
                 end
             end

             % check convergence
             if ~converged

                 % per-algorithm convergence
                 if alg_conv

                     % update counter
                     if abs(output.previous_best_funcvalues(stream,algorithm) - ...
                             output.best_funcvalues(stream,algorithm)) <= options.TolFun &&...
                        all(abs(output.previous_best_individuals(stream,algorithm,:) - ...
                             output.best_individuals(stream,algorithm,:)) <= options.TolX)
                         counter = counter + 1;
                     else counter = 0;
                     end

                     % if in addition, the minimum function value is less than
                     % options.QuitWhenAchieved is selected, we have convergence
                     if options.QuitWhenAchieved && ...
                        ((~constrained && output.global_best_funval <= options.AchieveFunVal) ||...
                         ( constrained && output.global_best_funval <= options.AchieveFunVal &&...
                           abs(output.global_best_funval_constrviolation) <= options.TolCon))
                         converged = true;
                         output.exitflag = 1;
                         output.message = sprintf(['Optimization terminated:\n\n',...
                             ' Best function value is less than OPTIONS.AchieveFunVal, and\n',...
                             ' OPTIONS.QuitWhenAchieve is switched ON.']);
                         return
                     end

                     % if counter is larger than preset maximum,
                     % convergence has been achieved
                     if (counter > options.TolIters)
                         converged = true;
                         % also set the output message and exitflag
                         % when there is only one stream
                         if number_of_streams == 1
                             output.exitflag = 1;
                             output.message =  sprintf(['Optimization terminated:\n\n',...
                             ' Coordinate differences were less than OPTIONS.TolX, and decrease\n',...
                             ' in function value was less than OPTIONS.TolFun for %d consecutive\n',...
                             ' iterations of the %s algorithm.\n',...
                             ' Note that since only 1 stream was used, further optimization is \n',...
                             ' impossble, so the demand of OPTIONS.MinIters of %d could not be met.\n'],...
                             options.TolIters, pop{i,j}.algorithm, options.MinIters);
                         end
                         % and return
                         return
                     end

                     % also, the standard deviation of function values in a
                     % population may decrease below options.TolX. In those
                     % cases, it's more efficient to just give up and
                     % continue with the next stream/algorithm
                     if std(pop{stream,algorithm}.fitnesses) < options.TolFun
                         converged = true;
                         % also set the output message and exitflag
                         % when there is only one stream
                         if number_of_streams == 1
                             output.exitflag = 1;
                             output.message =  sprintf(['Optimization terminated:\n\n',...
                             ' Standard deviation of function values was less than OPTIONS.TolFun.\n',...
                             ' Since only 1 stream was used, further optimization is not likely to\n',...
                             ' improve solutions, so GODLIKE decided to quit.\n']);
                         end
                         % and return
                         return;
                     end

                 % GODLIKE-convergence
                 else
                     % update counter
                     if output.global_best_funval < options.AchieveFunVal
                         if abs(output.previous_global_best_funval - ...
                                 output.global_best_funval) <= options.TolFun && ...
                            all(abs(output.previous_global_best_individual - ...
                                 output.global_best_individual) <= options.TolX)
                             output.descent_counter = output.descent_counter + 1;
                         else output.descent_counter = 0;
                         end
                     end % if

                     % if counter is larger than preset maximum, and the
                     % minimum amount of iterations has been performed,
                     % convergence has been achieved
                     if generation > options.MinIters && (output.descent_counter > 2)
                         converged = true;
                         output.exitflag = 1;
                         output.message = sprintf(['Optimization terminated:\n\n',...
                         ' Coordinate differences were less than OPTIONS.TolX, and decrease\n',...
                         ' in function value was less than OPTIONS.TolFun for 2 consecutive\n',...
                         ' GODLIKE-iterations. A total of %d iterations have been performed, \n',...
                         ' which satisfies OPTIONS.Miniters setting of %d.\n'],...
                         generation, options.MinIters);
                     end % if
                 end % if

                 % finalize output
                 if converged && ~alg_conv
                     % insert the last population in the output
                     for ii = 1:number_of_streams
                         % find the last non-empty population in this stream
                         jj = sum(~cellfun(@isempty, pop(ii, :)));
                         % adjust the output
                         if number_of_streams > 1
                             output.algorithm_info.(...
                                 ['Stream',num2str(ii)]).(upper(which_ones{jj}))(...
                                 sames(jj)).last_population = pop{ii,jj}.individuals;
                             output.algorithm_info.(...
                                 ['Stream',num2str(ii)]).(upper(which_ones{jj}))(...
                                 sames(jj)).last_fitnesses = pop{ii,jj}.fitnesses;
                         else
                             output.algorithm_info.(upper(which_ones{jj}))(...
                                 sames(jj)).last_population = pop{ii,jj}.individuals;
                             output.algorithm_info.(upper(which_ones{jj}))(...
                                 sames(jj)).last_fitnesses = pop{ii,jj}.fitnesses;
                         end
                     end
                 end
             end % if
         end % if single

         % stuff specific for multi-objective optimization
         if multi

             % check convergence
             if ~converged
                 % check if THIS algorithm has converges
                 if alg_conv
                     % are all individuals non-dominated?
                     converged = all(pop{stream,algorithm}.pop_data.front_number == 0); return

                 % check if ALL populations are non-dominated
                 else
                     all_nd = false(number_of_streams,1);
                     for ii = 1:number_of_streams
                         % find the last non-empty population in this stream
                         jj = sum(~cellfun(@isempty, pop(ii, :)));
                         % population might not have performed any iterations yet
                         if (pop{ii,jj}.iterations == 0), all_nd(ii) = true; continue; end
                         % otherwise, check the frontnumbers
                         all_nd(ii) = all(pop{ii,jj}.pop_data.front_number == 0);
                     end
                     converged = all(all_nd);
                 end
             end % if (~converged)

             % if converged, complete output
             if converged && ~alg_conv && ...
               (number_of_streams == 1 || generation > options.MinIters ||...
               generation >= options.MaxIters || num_funevaluations >= options.MaxFunEvals)

                 % Get Pareto-solutions & function values, and the "most" efficient one
                 [output.pareto_front_individuals, output.pareto_front_fitnesses, ...
                 output.most_efficient_point, output.most_efficient_fitnesses] = get_Paretos;

                 % and finalize the output structure
                 if (generation > options.MinIters) && (generation < options.MaxIters) &&...
                        (num_funevaluations <= options.MaxFunEvals)
                     output.exitflag = 1;
                     output.message = sprintf(['Optimization terminated:\n\n',...
                         ' All trial solutions of all selected algoritms are non-dominated.\n',...
                         ' A total of %d iterations have been performed, \n',...
                         ' which satisfies OPTIONS.MinIters setting of %d.\n'],...
                         generation, options.MinIters);
                 elseif (number_of_streams == 1) && (generation < options.MaxIters) &&...
                        (num_funevaluations <= options.MaxFunEvals)
                     % complete output structure
                     output.exitflag = 1;
                     output.message = sprintf(['Optimization terminated:\n',...
                         ' All members of the population are now non-dominated.\n',...
                         ' Note that since only 1 stream was used, further optimization is \n',...
                         ' impossble, so the demand of OPTIONS.MinIters of %d could not be met.'],...
                         options.MinIters);
                 end % if (MinIters check)

             % the minimum amount of iterations MUST be performed
             else
                 converged = false;

             end % if converged

         end % if multi

    end % test convergence

    % form the proper [x] and [optimValues] for output functions
    function [x, optimValues] = get_outputFcn_values(stream, algorithm)
        % collect all the information one could possibly desire
        if ~isempty(stream) && ~isempty(algorithm)
            optimValues.optimizer.algorithm  = pop{stream,algorithm}.algorithm;
            optimValues.optimizer.funcCount  = pop{stream,algorithm}.funevals;
            optimValues.optimizer.iterations = pop{stream,algorithm}.iterations;
            optimValues.stream = stream;
        else
            optimValues.optimizer.algorithm  = [];
            optimValues.optimizer.funcCount  = [];
            optimValues.optimizer.iterations = [];
            optimValues.stream = number_of_streams;
        end
        optimValues.algorithm = 'GODLIKE';
        optimValues.funcCount = num_funevaluations;
        optimValues.iteration = generation;
        optimValues.popsize   = popsize;
        if single
            optimValues.type = 'single-objective';
            if constrained
                optimValues.best_fval = output.global_best_funval_unconstrained;
                optimValues.best_fval_constraint_violation = ...
                    output.global_best_funval_constrviolation;
            else
                optimValues.best_fval = output.global_best_funval;
            end
            x = reshape(output.global_best_individual, sze);
            optimValues.best_individual = x;
        else
            optimValues.type = 'multi-objective';
            [non_dominated, non_dominated_fits, ...
                most_efficient_point, most_efficient_fitnesses, violations] = get_Paretos();
            optimValues.number_of_nondominated_solutions = size(non_dominated,1);
            optimValues.nondominated_solutions = non_dominated;
            optimValues.nondominated_function_values = non_dominated_fits;
            optimValues.most_efficient_point = reshape(most_efficient_point, sze);
            optimValues.most_efficient_fitnesses = most_efficient_fitnesses;
            x = reshape(most_efficient_point, sze);
            if constrained
                optimValues.constraint_violations = violations;
            end
        end
    end

    % Get the complete Pareto-front, and the most efficient point
    % (multi-objective optimization only)
    function [non_dominated, non_dominated_fits, ...
              most_efficient_point, most_efficient_fitnesses, violations] = get_Paretos
        % for unconstrained problems
        violations = [];
        % get complete current Pareto front
        non_dominated = [];
        non_dominated_fits = [];
        for s = 1:number_of_streams
            % find the last non-empty population in this stream
            t = sum(~cellfun(@isempty, pop(s, :)));
            % the algorithm might not have been used yet
            if (pop{s,t}.iterations == 0), continue; end
            % get the indices for the non-dominated members
            Pareto_indices = (pop{s,t}.pop_data.front_number == 0);
            % append its individuals
            non_dominated = ...
                [non_dominated
                pop{s,t}.individuals(Pareto_indices, :)];%#ok
            % and its fitnesses
            if constrained
                non_dominated_fits = ...
                    [non_dominated_fits
                     pop{s,t}.pop_data.unpenalized_function_values_parent(Pareto_indices, :)];%#ok
                violations = ...
                    [violations
                     pop{s,t}.pop_data.constraint_violations_parent(Pareto_indices, :)];%#ok
            else
                non_dominated_fits = ...
                    [non_dominated_fits
                    pop{s,t}.fitnesses(Pareto_indices, :)];%#ok
            end
        end
        % find most efficient point and its fitnesses
        origin                   = min(non_dominated_fits);
        shifted_fitnesses        = bsxfun(@minus, non_dominated_fits, origin);
        ranges                   = max(non_dominated_fits) - origin;
        scaled_fitnesses         = bsxfun(@rdivide, shifted_fitnesses*min(ranges), ranges);
        distances_sq             = sum(scaled_fitnesses.^2,2);
        [mindist_sq, index]      = min(distances_sq);%#ok
        most_efficient_point     = non_dominated(index, :);
        most_efficient_fitnesses = non_dominated_fits(index, :);
    end % get Paretos

    % display the algorithm's progress
    function display_progress

        % current loop indices
        algorithm_iteration = k; algorithm = j; stream = i;

        % 'iter' is the only display option for now
        if strcmpi(options.display, 'iter')

            % if not converged, display all relevant information
            % from the current iteration
            if ~converged

                % display header if this is the first generation, first
                % algorithm and first algorithm iteration
                if (generation == 1) && (algorithm_iteration == 1) && (algorithm == 1) && (stream == 1)

                    % display which number_of_algorithms
                    if number_of_algorithms == 1
                        strings = '%s population.\n';
                    elseif number_of_algorithms == 2
                        strings = '%s and %s populations.\n';
                    else
                        strings = repmat('%s, ', 1, number_of_algorithms-1);
                        % insert newlines if neccessary
                        if number_of_algorithms > 4
                            for ii = 16:64:length(strings)
                                strings(ii:ii+4) = '\n%s,';
                            end
                        end
                        strings = [strings, 'and %s populations.\n'];
                    end

                    % output header
                    fprintf(1, ['\nGODLIKE optimization started with ', strings], which_ones{:});

                    % display single or multi-objective optimization, and
                    % population size, iterations low and high
                    if     single
                        fprintf(1,...
                           ['Performing single-objective optimization, with total population\n'...
                            'size of %d individuals. Lower bounds on algorithm iterations\n', ...
                            'is %d, upper bound is %d.\n'], popsize, options.ItersLb, ...
                            options.ItersUb);
                    elseif multi
                        fprintf(1,...
                           ['Performing multi-objective optimization, with %d objectives.\n',...
                            'Total population size is %d individuals. Lower bounds on\nalgorithm ',...
                            'iterations is %d, upper bound is %d.\n'], options.num_objectives,...
                            popsize, options.ItersLb, options.ItersUb);
                    end % if
                end % if

                % display new iteration
                if (algorithm == 1) && (stream == 1) && (algorithm_iteration == 1)
                    fprintf(1,...
                        ['\n###############################################################\n',...
                        '                         ITERATION %d\n'],generation);
                    if single
                        fprintf(1, ...
                            '          Current global minimum: %+1.8e\n',...
                            output.global_best_funval);
                    end
                    fprintf(1, ...
                        '###############################################################\n');

                end % if

                % display new stream
                if (algorithm == 1) && (algorithm_iteration == 1) && (number_of_streams > 1)
                    fprintf(1,...
                        ['\n===============================================================\n',...
                         '                         Stream %d                                 ',...
                         '\n===============================================================\n'],stream);
                end

                % display new algorithm's header
                if (algorithm_iteration == 1)
                    fprintf(1,...
                        ['· · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·\n',...
                        '                     %s algorithm, %s pass\n',...
                        '               popsize: %d, max.iterations: %d\n'],...
                        which_ones{algorithm}, counter_appendix(generation), ...
                        frac_popsize(stream), frac_iterations(algorithm));
                    if ~constrained
                        if multi
                            fprintf(1, ...
                                '  #  f.count      Pareto fronts       non-Pareto fronts\n');
                        elseif single
                            fprintf(1, '  #  f.count       min(F)        std(F)         descent\n');
                        end % if
                        fprintf(1, '· · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·\n');
                    else
                        if multi
                            fprintf(1, ...
                                '  #  f.count      Pareto fronts       non-Pareto fronts\n');
                        elseif single
                            fprintf(1, '  #  f.count       min(F)        min(F_con)   min(pen.(F))      descent\n');
                        end % if
                        fprintf(1, '· · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·\n');
                    end
                end
                if ~constrained
                    if multi
                        fprintf(1, '%3d   %6d    %10d             %10d\n', ...
                            algorithm_iteration, pop{stream,algorithm}.funevals, ...
                            nnz(pop{stream,algorithm}.pop_data.front_number==0),...
                            nnz(pop{stream,algorithm}.pop_data.front_number~=0));
                    elseif single
                        fprintf(1, '%3d   %6d    %+1.5e  %+1.5e  %+1.5e\n',...
                            algorithm_iteration, pop{stream,algorithm}.funevals, ...
                            min(pop{stream,algorithm}.fitnesses),...
                            std(pop{stream,algorithm}.fitnesses),...
                            output.previous_best_funcvalues(stream,algorithm) -...
                            output.best_funcvalues(stream,algorithm));
                    end % if
                else
                    if multi
                        fprintf(1, '%3d   %6d    %10d             %10d\n', ...
                            algorithm_iteration, pop{stream,algorithm}.funevals, ...
                            nnz(pop{stream,algorithm}.pop_data.front_number==0),...
                            nnz(pop{stream,algorithm}.pop_data.front_number~=0));
                    elseif single
                        [best_penalized_F, ind] = min(pop{stream,algorithm}.fitnesses);
                        best_unpenalized_F = ...
                            pop{stream,algorithm}.pop_data.unpenalized_function_values_parent(ind);
                        best_constraint_violation = ...
                            pop{stream,algorithm}.pop_data.constraint_violations_parent(ind);
                        fprintf(1, '%3d   %6d    %+1.5e  %+1.5e  %+1.5e  %+1.5e\n',...
                            algorithm_iteration, ...
                            pop{stream,algorithm}.funevals, ...
                            best_unpenalized_F,...
                            best_constraint_violation,...
                            best_penalized_F,...
                            output.previous_best_funcvalues(stream,algorithm) -...
                            output.best_funcvalues(stream,algorithm));
                    end % if
                end

            % if we do have convergence, just display the output message
            else
                fprintf(1, '\n'); fprintf(1, output.message); fprintf(1, '\n\n');
            end % if

        end % ITER

        % create proper counter and appendix
        % (as in 2 -> '2nd',  3 -> '3rd', etc.)
        function out_string = counter_appendix(integer)
            % create counter string
            inp_string = num2str(integer);
            if strcmp(inp_string,'11')||strcmp(inp_string,'12')||strcmp(inp_string,'13')
                out_string = [inp_string,'th'];
            else
                switch inp_string(end)
                    case '1', out_string  = [inp_string,'st'];
                    case '2', out_string  = [inp_string,'nd'];
                    case '3', out_string  = [inp_string,'rd'];
                    otherwise, out_string = [inp_string,'th'];
                end
            end
        end

    end % display progress

%     % coordinate transformation
%     function trans_individuals = transform(individuals)
%
%         % get the indices for the fixed variables
%         eq_indices = pop{1}.eq_indices(1, :);
%
%         % remove fixed variables
%         trans_individuals = individuals(:, ~eq_indices);
%
%         % apply asin-transformation
%         trans_individuals = real(asin(2*(trans_individuals-lb(~eq_indices))./...
%             (ub(~eq_indices)-lb(~eq_indices))-1));
%     end % transform

end % GODLIKE
