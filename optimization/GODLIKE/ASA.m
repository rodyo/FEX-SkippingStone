classdef ASA < pop_multi
% DIFFERENTIAL EVOLUTION
%
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

    % public methods
    methods (Access = public)
        
        % simple constructor: create pop_single object, and 
        % just add the number of objectives
        function pop = ASA(varargin)
            % construct pop_single object
            pop = pop@pop_multi(varargin{:});            
            % ASA also requires additional initializations
            pop.initialize_algorithm;
            
        end % constructor
        
        % single iteration
        function iterate(pop)
            
            % the sequence of one iteration depends on whether the
            % optimization is single or multi-objective
            if pop.num_objectives == 1 % single-objective
                pop.create_offspring(1:pop.size);    % create offspring                
                pop.evaluate_function;               % evaluate objective function
                pop.replace_parents;                 % replace the parents                                
                
            else % multi-objective
                pop.non_dominated_sort;              % non-dominated sort
                pop.update_algorithm;                % update temperature
                pool = ...                           % binary tournament selection (full pop.size)
                    pop.tournament_selection_multi(pop.size, 2);
                pop.create_offspring(pool);          % create new offspring
                pop.evaluate_function;               % carefully evaluate objective function                
                
            end
            
             % increase number of iterations made
            pop.iterations = pop.iterations + 1;
                        
        end % single iteration
        
    end % public methods
    
    % hidden methods
    methods (Hidden)
        
        % generate new generation
        function create_offspring(pop, pool)
                        
            % rename some stuff 
            parent_pop = pop.trans_individuals(pool, :);
            parent_fit = pop.fitnesses(pool, :);
            reinit = pop.options.ReinitRatio;
            
            % initialize placeholders
            newpop = NaN(pop.size, pop.dimensions);
            newfit = NaN(pop.size, pop.options.num_objectives);     

            % first, make sure the pool is large enough
            temp_pop = zeros(pop.size, pop.dimensions);
            temp_fit = temp_pop;
            if numel(parent_pop) ~= pop.size*pop.dimensions
                % insert all values from the pool
                temp_pop(1:size(parent_pop,1), :) = parent_pop; 
                temp_fit(1:size(parent_pop,1), :) = parent_fit; 
                % insert random members from [parent_pop]
                for i = size(parent_pop,1)+1:pop.size 
                    randinds = round(rand*(size(parent_pop,1)-1)+1); 
                    temp_pop(i, :) = parent_pop(randinds, :);
                    temp_fit(i, :) = parent_fit(randinds, :);
                end
                % equate the two
                parent_pop = temp_pop;
                parent_fit = temp_fit;
            end
            
            % reinitialize a small percentage of the population  
            reinit = floor(pop.size * reinit);
            if (reinit > 0)
                % reinitialize
                new_individuals = pop.lb(1:reinit,:) + ...
                        rand(reinit, pop.dimensions) .* (pop.ub(1:reinit,:)-pop.lb(1:reinit,:));
                % apply sine-transformation
                new_individuals = real(asin(2*(new_individuals-pop.lb(1:reinit,:))./...
                    (pop.ub(1:reinit,:)-pop.lb(1:reinit,:))-1));
                % Use Bolzmann distribution to create new points
                rands = sqrt(pop.pop_data.temperature)*randn(pop.size-reinit, pop.dimensions);
                
                % Use worst percentage for single-objective 
                if (pop.dimensions == 1)
                    % sort population
                    [sorted_fits, indices] = sort(parent_fit, 1, 'descend');%#ok
                    % reinitialize
                    newpop(indices(1:reinit),:) = new_individuals;
                    newfit(indices(1:reinit),:) = NaN;                    
                    newpop(indices(reinit+1:end),:) = parent_pop(indices(reinit+1:end),:) + rands;
                % use random percentage for multi-objective
                else
                    % reinitialize
                    newpop(1:reinit,:) = new_individuals;
                    newfit(1:reinit,:) = NaN;                    
                    newpop(reinit+1:end,:) = parent_pop(reinit+1:end,:) + rands;
                end
            end            
            
            % at least one dimension always changes, so ALL function
            % values have to be recomputed
            
            % insert result into pop
            pop.pop_data.offspring_population      = newpop;
            pop.pop_data.function_values_offspring = newfit;            
            
        end
        
        % selectively replace the parent population with members
        % from the offspring population (single-objective optimization)
        function replace_parents(pop)  
            
            % rename some stuff                    
            new_fits = pop.pop_data.function_values_offspring;
            new_inds = pop.pop_data.offspring_population;
            T        = pop.pop_data.temperature;
            T0       = pop.options.ASA.T0;
            nrg      = pop.pop_data.function_values_offspring;
            prevnrg  = pop.pop_data.function_values_parent;
            cool     = pop.options.ASA.CoolingSchedule;
            iters    = pop.iterations - pop.pop_data.iters;
            
            % reject or accept the new population, according to
            % the probabilistic rule
            nrgdiff  = (prevnrg - nrg);                % energy difference
            nrgdiff  = nrgdiff / max(abs(nrgdiff(:))); % rescale the differences
            ind      = nrgdiff > 0;                    % always accept better ones
            probind  = ~ind & rand(pop.size, 1) < exp( nrgdiff/T );
            % accept worse ones based on
            % probabalistic rule
            swapinds = ind|probind;                    % indices to be swapped
            
            % apply cooling schedule
            pop.pop_data.temperature = max(eps,cool(T, T0, iters));
            
            % replace the individuals & function values
            pop.pop_data.parent_population(swapinds, :)      = new_inds(swapinds, :);
            pop.pop_data.function_values_parent(swapinds, :) = new_fits(swapinds, :);
            % constraint functions
            if ~isempty(pop.constrained)
                pop.pop_data.unpenalized_function_values_parent(swapinds, :) = ...
                    pop.pop_data.unpenalized_function_values_offspring(swapinds, :);
                pop.pop_data.constraint_violations_parent(swapinds, :) = ...
                    pop.pop_data.constraint_violations_offspring(swapinds, :);
            end
            
            % copy individuals and fitnesses to respective properties            
            pop.fitnesses = pop.pop_data.function_values_parent;       
            pop.trans_individuals = pop.pop_data.parent_population;   
            pop.individuals = pop.wrapperFcn(pop.trans_individuals, 1:pop.size);
            
        end % replace parents
        
        % initialize the algorithm
        function initialize_algorithm(pop)
            
            % if the initial temperature is left empty, estimate
            % an optimal one. This is simply the largest quantity
            % in (ub - lb), divided by 4, and squared. This
            % ensures that during the first few iterations,
            % particles are able to spread over the entire search
            % space; 4 = 2*2*std(randn(inf,1)).
            if isempty(pop.options.ASA.T0)
                % only do it upon initialization
                
                % find the maximum value
                sqrtT0dv4 = max(pop.ub(1,:) - pop.lb(1, :))/4;
                
                % set T0
                pop.options.ASA.T0 = sqrtT0dv4^2;
                
            end
            
            % initialize temperature
            pop.pop_data.temperature = pop.options.ASA.T0;
            
            % initialize iterations
            pop.pop_data.iters = pop.iterations;
            
            
        end % initialize algorithm
        
        % Update globally changing variables associated with PSO
        function update_algorithm(pop)
            
            % rename some stuff
            T       = pop.pop_data.temperature;
            T0      = pop.options.ASA.T0;
            cool    = pop.options.ASA.CoolingSchedule;
            iters   = pop.iterations - pop.pop_data.iters;
            
            % essentially, only the temperature is lowered
            pop.pop_data.temperature = max(eps, cool(T, T0, iters));
            
        end % update algorithm
        
    end % protected/hidden methods
    
end % classdef