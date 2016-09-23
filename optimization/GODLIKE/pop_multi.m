classdef pop_multi < pop_single
% POP_MULTI         Class definition for a population to be
%                   used for multi-objective optimization 
%
% POP_MULTI is a SubClass of POP_SINGLE. The class 
% constructor works in the same way as that of POP_SINGLE,
% with the exception that an additional property is set:
%
%   pop.num_objectives      (number of objectives)
%
% All inputs and other properties are the same as for 
% POP_SINGLE -- type 'help pop_single' for more information. 
%
% The method ITERATE is now suited to optimize multi-
% objective problems. To that end, several other (hidden) 
% methods have been implemented: 
%
%   NON_DOMINATED_SORT()
%
%       a general implementation of NSGA-II's non-dominated 
%       sorting procedure. Sorts the current population
%       according to the domination level of the individuals.
%
%
%   pool = TOURNAMENT_SELECTION(pool_size, tournament_size) 
%
%       a general tournament selection procedure, that takes
%       [tournament_size] individuals randomly selected from
%       the offspring population and lets them compete with
%       the rankings and crowding distances as determining 
%       factors. The winning individual of each tournament
%       is inserted into [pool] until that [pool] contains
%       [pool_size] individuals. 
%
%   UPDATE_ALGORITHMS()
%    
%       Called from NON_DOMINATED_SORT(), updates some 
%       globaly changing values for the different algorithms. 
%       In pop_single, this is done in REPLACE_PARENTS(), but
%       as that step is not executed here, an extra method is
%       required. This updates for instance the [lbest],
%       [nbest] and [gbest] for PSO, and the temperature for
%       ASA. 
%
%
% See also pop_single, GODLIKE. 

%      Author : Rody P.S. Oldenhuis
% Affiliation : Delft University of Technology
%               Faculty of Aerospace Engineering
%               Dep. of Astrodynamics & Satellite Systems 
%     Contact : oldnhuis@dds.nl
%   Licensing/
%    (C) info : Frankly I don't care what you do with it, 
%               as long as I get some credit when you copy 
%               large portions of the code ^_^

% last edited 03/10/2009
    
    % additional properties 
    properties
        num_objectives     % number of objectives
        % contents of pop_data for multi-objective
        % optimization:
        %      pop_data.parent_population
        %      pop_data.offspring_population
        %      pop_data.function_values_parent
        %      pop_data.function_values_offspring
        %      pop_data.front_number
        %      pop_data.crowding_distance        
    end
    
    % public methods
    methods (Access = public)
                
        % simple constructor: create pop_single object, and 
        % just add the number of objectives
        function pop = pop_multi(varargin)
            % construct pop_single object
            pop = pop@pop_single(varargin{:});
            % just set the amount of objectives
            pop.num_objectives = pop.options.num_objectives;
        end % constructor
        
    end % public methods
    
    % protected/hidden methods
    methods (Access = protected, Hidden)
        
        % non-dominated sort, and crowding distance assignment
        function non_dominated_sort(pop)
            
            % this function sorts the population according to non-domination.
            % At the same time, it will compute the crowding distance for every
            % individual. It will store these values for all individuals in the 
            % data structure [pop.pop_data]. 
                        
            % determine if this is the first call or a subsequent call, 
            % and create arrays accordingly
            if ~isempty(pop.pop_data.function_values_offspring)
                inds = [pop.pop_data.parent_population
                    pop.pop_data.offspring_population; ];
                fits = [pop.pop_data.function_values_parent
                    pop.pop_data.function_values_offspring];
                if ~isempty(pop.constrained)
                    cons = [pop.pop_data.constraint_violations_parent;
                        pop.pop_data.constraint_violations_offspring];
                end
                N = 2*pop.size;
            else
                inds = pop.trans_individuals;
                fits = pop.fitnesses;
                if ~isempty(pop.constrained)
                    cons = pop.constraint_violations_offspring;
                end
                N = pop.size;
            end
            
            % copy fitnesses            
            fits_copy = fits;
                             
            % pre-calculate some stuff for crowding distances
            crowding_dists = zeros(N, 1);              % initialize            
            [sorted_fitnesses, indices] = sort(fits);  % sort fitnesses            
            crowding_dists(indices(1, :))   = inf;     % always include boundaries
            crowding_dists(indices(end, :)) = inf;     % always include boundaries
              
            % initialize some variables
            guy_is_dominated_by = (2*pop.size+1)*ones(N, 1);
            front_number = 0;            
            
            % compute the front numbers of all individuals      
            % NOTE: two different routines: One for constrained
            % optimization, and one for unconstrained optimization.
            %
            % paretofront MEX-file by Yi Cao
            % (see MATLAB File Exchange, file 17251)
            while ~all(isinf(fits_copy))  
                %do the Pareto membersip test
                front = paretofront(fits_copy);
                %fill the array
                guy_is_dominated_by(front) = front_number;
                %make sure these individuals don't dominate anything in
                %the next iteration
                fits_copy(front, :) = inf;
                %increase the front number
                front_number = front_number + 1;
            end
                        
            % compute the crowding distances
            for guy = 1:N
                % compute this guy's crowding distance
                for m = 1:pop.num_objectives
                    % current sorting indices
                    sort_inds = indices(:, m);                    
                    % find this guy's index
                    guys_ind = find(sort_inds == guy);                    
                    % compute crowding distance
                    if isfinite(crowding_dists(guy))
                        crowding_dists(guy) = crowding_dists(guy) + ...
                            (sorted_fitnesses(guys_ind+1, m)-sorted_fitnesses(guys_ind-1, m))/...
                            (sorted_fitnesses(end,m)-sorted_fitnesses(1,m));
                    else break
                    end % if
                end % for                
            end % for
            
            % create new population
            new_pop = zeros(pop.size, pop.dimensions);
            new_fit = zeros(pop.size, pop.num_objectives);
            fronts  = zeros(pop.size, 1);
            not_filled = pop.size;
            % also handle constraints
            if ~isempty(pop.constrained)
                new_cons = zeros(pop.size, size(cons, 2));
            end
            for i = 0:max(guy_is_dominated_by)
                % extract current front
                front = guy_is_dominated_by == i;                
                % number of entries
                entries = nnz(front);                
                % see if it still fits in the population 
                if entries <= not_filled                     
                    % if so, insert all individuals from this front in the population
                    % and keep track of their fitnesses                    
                    new_pop(pop.size-not_filled+1:pop.size-not_filled+entries, :) = inds(front, :);
                    new_fit(pop.size-not_filled+1:pop.size-not_filled+entries, :) = fits(front, :);
                    fronts (pop.size-not_filled+1:pop.size-not_filled+entries, 1) = i;
                    % handle constraints
                    if ~isempty(pop.constrained)
                        new_cons(pop.size-not_filled+1:pop.size-not_filled+entries, :) = cons(front, :);
                    end
                    % adjust number of entries that have not yet been filled
                    not_filled = not_filled - entries;                  
                % if it does not fit, insert the remaining individuals based on 
                % their crowding distance
                else break 
                end % if                
            end % for
              
            % (for the first iteration, the WHOLE current population will fit 
            % in the new population. Skip that case)
            if (N ~= pop.size)
                % sort last front encountered by crowding-distance
                front_inds = inds(front, :);       front_fits = fits(front, :);
                [sorted_front, indices] = sort(crowding_dists(front), 'descend');%#ok             
                % extract individuals & fitnesses in proper order
                sorted_inds = front_inds(indices, :); sorted_fits = front_fits(indices, :);                    
                % insert the remaining individuals in the new population 
                new_pop(pop.size-not_filled+1:end, :) = sorted_inds(1:not_filled, :);
                new_fit(pop.size-not_filled+1:end, :) = sorted_fits(1:not_filled, :);
                fronts (pop.size-not_filled+1:end, 1) = i; 
            end % if 

            % insert results in data structure                
            pop.pop_data.parent_population      = new_pop;            
            pop.pop_data.function_values_parent = new_fit;
            pop.pop_data.front_number           = fronts;
            pop.pop_data.crowding_distance      = crowding_dists;  
            
            % copy individuals and fitnesses to respective class properties
            pop.fitnesses = new_fit;
            pop.trans_individuals = new_pop;
            pop.individuals = pop.wrapperFcn(new_pop, 1:pop.size);
            % also copy fitnesses to appropriate field for constrained
            % optimizations
            if ~isempty(pop.constrained)
                pop.pop_data.unpenalized_function_values_parent = pop.fitnesses;
            end
            
        end % non-dominated sort
        
        % tournament selection with crowding distances and rankings 
        % as competitive factors
        function pool = tournament_selection_multi(pop, pool_size, tournament_size)
            
            % initialize mating pool
            pool = zeros(pool_size, 1);
            
            % fill the mating pool
            for i = 1:pool_size    
                
                % select [tournament_size] individuals at random
                equal = true;
                while equal % make sure the indices are not equal
                    inds  = round( rand(tournament_size, 1)*(pop.size-1) + 1); % random indices
                    equal = numel(unique(inds)) ~= numel(inds);        % check if any are equal
                end
                
                % let them compete according to
                % (xj < yj) if (rank(xj) < rank(yj))
                %           or (rank(xj) == rank(yj) and distance(xj) > distance(yj)
                %
                % The definition of [rank] is different for unconstrained and
                % constrained problems:
                %
                % UNCONSTRAINED PROBLEMS                
                % [rank] is simply the front number of that solution
                %
                % CONSTRAINED PROBLEMS                
                %    rank(xj) < rank(yj) if 
                %       1) Solution (xj) is feasible and solution (yj) is not
                %       2) Both solutions are infeasible, but (xj) has a
                %          smaller overall constraint violation than (yj)
                %       3) Both solutions are feasible, and (xj) dominates
                %          solution (yj)
                ranks     = pop.pop_data.front_number(inds);
                distances = pop.pop_data.crowding_distance(inds);
                if ~isempty(pop.constrained)
                    con_violations = pop.pop_data.constraint_violations_parent(inds);
                end
                for j = 1:tournament_size 
                    % find those with lower constraint violations
                    if ~isempty(pop.constrained)
                        % compare constraint violations
                        less_con_violation = con_violations(j) < ...
                            [con_violations(1:j-1); con_violations(j+1:end)];
                        % constr. violation is less than all others
                        % (implements (1) and (2) simultaneously)
                        if all(less_con_violation), best = inds(j); break; end
                    end
                    % find those of lower rank
                    % (implements (3) of constrained optimization)
                    less_rank = ranks(j) < [ranks(1:j-1); ranks(j+1:end)];
                    % rank is less than all others
                    if all(less_rank), best = inds(j); break; end
                    % compare distances and rank equality
                    more_dist = ranks(j) == [ranks(1:j-1); ranks(j+1:end)] &...
                                distances(j) >= [distances(1:j-1); distances(j+1:end)];                    
                    % rank is equal, but distance is less than all others
                    if any(~less_rank & more_dist), best = inds(j); break; end      
                    % if no best has been found yet, the solutions are all
                    % equal. Just output the current solution in that case
                    best = inds(j);                    
                end % for
                
                % insert the index of the best one in the pool
                pool(i) = best;                
            
            end % for
            
        end % tournament selection multi
                        
    end % protected/hidden methods 
    
end % classdef
