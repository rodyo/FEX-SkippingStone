classdef DE < pop_multi
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
        function pop = DE(varargin)
            % construct pop-object
            pop = pop@pop_multi(varargin{:});
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
                pool = ...                           % binary tournament selection (full pop.size)
                    pop.tournament_selection_multi(pop.size, 2);                
                pop.create_offspring(pool);          % create new offspring                
                pop.evaluate_function;               % evaluate objective function                                
            end
            
            % increase number of iterations made
            pop.iterations = pop.iterations + 1; 
                        
        end % single iteration              
        
    end % public methods
    
    % hidden methods
    methods (Hidden)
        
        % generate new generation
        function create_offspring(pop, pool)
            
            % get the size of the pool
            pool_size = length(pool);
            
            % rename some stuff 
            parent_pop = pop.trans_individuals(pool, :);
            parent_fit = pop.fitnesses(pool, :);
            
            % initialize
            newpop = zeros(pop.size, pop.dimensions);           % empty new population            
            newfit = NaN(pop.size, pop.options.num_objectives); % placeholder for the sites to 
                                                                % evaluate the function       
            % rename some stuff
            Flb = pop.options.DE.Flb;
            Fub = pop.options.DE.Fub;
            crossconst = pop.options.DE.CrossConst;
            reinit = pop.options.ReinitRatio;
            
            % reinitialize a small percentage of the population  
            reinit = floor(pop.size * reinit);
            if (reinit > 0)
                % reinitialize
                new_individuals = pop.lb(1:reinit,:) + ...
                        rand(reinit, pop.dimensions) .* (pop.ub(1:reinit,:)-pop.lb(1:reinit,:));
                % apply sine-transformation
                new_individuals = real(asin(2*(new_individuals-pop.lb(1:reinit,:))./...
                    (pop.ub(1:reinit,:)-pop.lb(1:reinit,:))-1));
                
                % Use worst percentage for single-objective 
                if (pop.dimensions == 1)
                    % sort population
                    [sorted_fits, indices] = sort(parent_fit, 1, 'descend');%#ok
                    % reinitialize
                    newpop(indices(1:reinit),:) = new_individuals;
                    newfit(indices(1:reinit),:) = NaN;
                % use random percentage for multi-objective
                else
                    newpop(1:reinit,:) = new_individuals;
                    newfit(1:reinit,:) = NaN;
                end
            end           
                        
            % Neoteric Differential Evolution
            for i = (reinit+1):pop.size
                % random indices
                base = round(rand*(pool_size-1))+1; % RANDI is slower
                d1   = round(rand*(pool_size-1))+1;
                d2   = round(rand*(pool_size-1))+1;
                % d2 may not be equal to d1
                while (d1 == d2), d2 = round(rand*(pool_size-1))+1; end
                % DE operator
                if (rand < crossconst) || round(rand*(pool_size-1))+1 == i;
                    % DE operator when rnd < Cr
                    F = rand*(Fub-Flb) + Flb;
                    newpop(i, :) = parent_pop(base,:) +  ...
                        F*(parent_pop(d1,:) - parent_pop(d2,:));
                else
                    % insert random parent otherwise
                    rnd_ind      = round(rand*(pool_size-1))+1;
                    newpop(i, :) = parent_pop(rnd_ind, :);
                    newfit(i, :) = parent_fit(rnd_ind, :);                    
                end
            end % for
            
            % insert result into pop
            pop.pop_data.offspring_population      = newpop;
            pop.pop_data.function_values_offspring = newfit;
                         
        end
        
        % selectively replace the parent population with members
        % from the offspring population (single-objective optimization)
        function replace_parents(pop)  
            
            % rename for clarity            
            new_fits   = pop.pop_data.function_values_offspring;
            old_fits   = pop.fitnesses;
            new_inds   = pop.pop_data.offspring_population;               
            
            % DE uses simple greedy replacement            
            better_inds = new_fits < old_fits;
            pop.pop_data.parent_population(better_inds, :) = new_inds(better_inds, :);
            pop.pop_data.function_values_parent(better_inds, :) = new_fits(better_inds, :);
            % constraint functions
            if ~isempty(pop.constrained)
                pop.pop_data.unpenalized_function_values_parent(better_inds, :) = ...
                    pop.pop_data.unpenalized_function_values_offspring(better_inds, :);
                pop.pop_data.constraint_violations_parent(better_inds, :) = ...
                    pop.pop_data.constraint_violations_offspring(better_inds, :);
            end  
            
            % copy individuals and fitnesses to respective properties
            pop.fitnesses = pop.pop_data.function_values_parent;
            pop.trans_individuals = pop.pop_data.parent_population;
            pop.individuals = pop.wrapperFcn( pop.trans_individuals, 1:pop.size);                
            
        end % replace parents
        
    end % protected/hidden methods
    
end % classdef