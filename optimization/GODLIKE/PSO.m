classdef PSO < pop_multi
% PARTICLE SWARM OPTIMIZATION
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
        function pop = PSO(varargin)
            % construct pop_single object
            pop = pop@pop_multi(varargin{:}); 
            % PSO also requires additional initializations
            pop.initialize_algorithm;
        end % constructor
        
        % single iteration
        function iterate(pop)
                     
            % the sequence of one iteration depends on whether the
            % optimization is single or multi-objective
            if pop.num_objectives == 1 % single-objective
                pop.create_offspring(1:pop.size); % create offspring
                pop.evaluate_function;            % evaluate objective function
                pop.replace_parents;              % replace the parents
                
            else % multi-objective
                pop.non_dominated_sort;           % non-dominated sort   
                pool = ...                        % binary tournament selection (full pop.size)
                    pop.tournament_selection_multi(pop.size, 2);
                pop.create_offspring(pool);       % create new offspring
                pop.evaluate_function;            % carefully evaluate objective function
            end
            
            % update neighbors, speeds etc.
            pop.update_algorithm;                   
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
            parent_velocities = pop.pop_data.velocities(pool, :);
            reinit = pop.options.ReinitRatio;
                        
            % first, make sure the pool is large enough
            temp_pop = zeros(pop.size, pop.dimensions);
            temp_vel = temp_pop;
            temp_fit = temp_pop;
            if size(parent_pop,1) ~= pop.size
                % insert all values from the pool
                temp_pop(1:size(parent_pop,1), :) = parent_pop;
                temp_vel(1:size(parent_pop,1), :) = parent_velocities;
                temp_fit(1:size(parent_pop,1), :) = parent_fit;
                % insert random members from [parent_pop]
                for i = size(parent_pop,1)+(1:pop.size) 
                    randinds = round(rand*(size(parent_pop,1)-1)+1); 
                    temp_pop(i, :) = parent_pop(randinds, :);
                    temp_vel(i, :) = parent_velocities(randinds, :);
                    temp_fit(i, :) = parent_fit(randinds, :);
                end
                % equate the two
                parent_pop = temp_pop;
                parent_velocities = temp_vel;
                parent_fit = temp_fit;
            end
            
            % Creating offspring with PSO is pretty simple:
            newpop = parent_pop + parent_velocities;
            
            % "teleport" a small percentage of the population  
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
                % use random percentage for multi-objective
                else
                    newpop(1:reinit,:) = new_individuals;
                end
            end
            
            % Since the velocity can not be zero, each individual
            % changes during the creation of offspring, so the
            % function values of ALL new individuals have to be
            % re-calculated.
            
            % insert result into pop
            pop.pop_data.offspring_population      = newpop;
            pop.pop_data.function_values_offspring = NaN(pop.size, pop.options.num_objectives);
            
        end
        
        % selectively replace the parent population with members
        % from the offspring population (single-objective optimization)
        function replace_parents(pop)  
            
            % rename for clarity
            new_fits = pop.pop_data.function_values_offspring;            
            new_inds = pop.pop_data.offspring_population;
            
            % PSO simply replaces all parents
            pop.pop_data.parent_population = new_inds;
            pop.pop_data.function_values_parent = new_fits;
            % constraint functions
            if ~isempty(pop.constrained)
                pop.pop_data.unpenalized_function_values_parent = ...
                    pop.pop_data.unpenalized_function_values_offspring;
                pop.pop_data.constraint_violations_parent = ...
                    pop.pop_data.constraint_violations_offspring;
            end
                        
            % copy individuals and fitnesses to respective properties            
            pop.fitnesses         = new_fits;       
            pop.trans_individuals = new_inds;   
            pop.individuals       = pop.wrapperFcn(new_inds, 1:pop.size);
            
        end % replace parents
        
        % initialize algorithm
        function initialize_algorithm(pop)
            
            % initialize velocities
            % (average velocities about 20% of [[lb]-[ub]] interval)
            pop.pop_data.velocities = randn(pop.size, pop.dimensions) .* (pop.ub-pop.lb)/5;
            
            % compute speed limits            
            lower_speed_limit = (pop.ub(1, :)-pop.lb(1, :))/1e14;
            pop.pop_data.lower_speed_limit = sqrt(lower_speed_limit*lower_speed_limit.');
            upper_speed_limit = (pop.ub(1, :)-pop.lb(1, :))/10;
            pop.pop_data.upper_speed_limit = sqrt(upper_speed_limit*upper_speed_limit.');
            
            % rename for clarity
            NumNeighbors = pop.options.PSO.NumNeighbors;
            NetworkTopology = pop.options.PSO.NetworkTopology;
            
            % initialize neighbors
            switch lower(NetworkTopology)
                
                case 'star'
                    % star topology - in each star, there is one
                    % focal particle, to which the other members
                    % of the star are connected.
                    
                    % initialize
                    all_particles = (1:pop.size).';
                    num_stars     = floor(pop.size/NumNeighbors);
                    pop.pop_data.neighbors = zeros(pop.size, NumNeighbors-1);
                    
                    % initialize stars
                    if num_stars ~= 0
                        [dummy, focals] = sort(rand(pop.size,1));%#ok
                        focals = all_particles(focals(1:num_stars));
                        all_particles(focals) = [];
                        % select [NumNeighbors] random & unique neighbors
                        % for each focal particle
                        for i = 1:num_stars
                            % select new neighbors
                            [dummy, inds] = sort(rand(size(all_particles,1),1));%#ok
                            new_neighs = all_particles(inds(1:NumNeighbors-1));
                            % adjust array
                            for j = 1:NumNeighbors-1
                                all_particles(all_particles == new_neighs(j)) = [];
                            end
                            % assign new neighbors to focal particle
                            pop.pop_data.neighbors(focals(i), :) = new_neighs;
                            % assign focal particle to new neighbors
                            pop.pop_data.neighbors(new_neighs, 1) = focals(i);
                        end
                    else
                    end
                    
                    % population might be badly scaled for selected
                    % number of stars. Correct for this
                    % TODO - it works; those particles simply have no
                    % neighbors
                    
                case 'ring'
                    
                    % initialize
                    all_particles = (1:pop.size).';
                    num_rings     = floor(pop.size/NumNeighbors);
                    pop.pop_data.neighbors = zeros(pop.size, 2);
                    
                    % form the ring
                    for i = 1:num_rings
                        % randomly select [NumNeighbors] particles
                        [dummy, inds] = sort(rand(size(all_particles,1),1));%#ok
                        new_neighs = all_particles(inds(1:NumNeighbors));
                        % insert circularly shifted arrays
                        pop.pop_data.neighbors(new_neighs, :) = ...
                            [circshift(new_neighs,1), circshift(new_neighs,-1)];
                        % adjust array
                        for j = 1:NumNeighbors
                            all_particles(all_particles == new_neighs(j)) = [];
                        end
                    end
                    
                    % population might be badly scaled for selected
                    % number of rings. Correct for this!
                    % TODO - it works; those particles simply have no
                    % neighbors
                    
                case 'fully_connected'
                    % fully connected swarm - all particles have
                    % ALL other particles as neighbor
                    
                    % initialize
                    pop.pop_data.neighbors = zeros(pop.size, pop.size-1);
                    
                    % fill the neighbors
                    for i = 1:pop.size
                        pop.pop_data.neighbors(i, :) = [1:i-1, i+1:pop.size];
                    end
                    
            end % switch
            
            % find global best solution
            [global_best, index] = min(pop.fitnesses(:,1));%#ok
            pop.pop_data.global_best_ind = pop.trans_individuals(index, :);
            pop.pop_data.global_best_fit = pop.fitnesses(index, :);
            
            % initially, local best solutions are the function values themselves
            pop.pop_data.local_best_inds = pop.trans_individuals;
            pop.pop_data.local_best_fits = pop.fitnesses;
            
            % find the neighbor best
            pop.pop_data.neighbor_best_fits = zeros(pop.size, pop.options.num_objectives);
            pop.pop_data.neighbor_best_inds = zeros(pop.size, pop.dimensions);
            for i = 1:pop.size
                neighbors = pop.pop_data.neighbors(i, :);
                neighbors = neighbors(neighbors ~= 0);
                if isempty(neighbors), continue, end
                [neighbor_best, ind] = min(pop.fitnesses(neighbors));
                pop.pop_data.neighbor_best_fits(i, 1) = neighbor_best;
                pop.pop_data.neighbor_best_inds(i, :) = pop.trans_individuals(ind, :);
            end
            
        end % initialize algorithm
        
        % Update globally changing variables associated with PSO
        function update_algorithm(pop)
            
            % single objective
            if pop.num_objectives == 1
                
                % rename for clarity
                new_fits = pop.pop_data.function_values_offspring;
                new_inds = pop.pop_data.offspring_population;
                
                % update the neighbor bests
                % (this implementation is fast, but not intuitive)
                % add one NaN to the new_fits array
                new_fits  = [new_fits; NaN];
                % copy neighbors
                neighbors = pop.pop_data.neighbors;
                % let those that are zero refer to the NaN entry
                neighbors(neighbors == 0) = size(new_fits,1);
                % find the best ones
                [neighbor_best, ind] = min(new_fits(neighbors),[],2);
                % find those that are better
                better_neighbors = neighbor_best <  pop.pop_data.neighbor_best_fits;
                % no better ones might be found
                if any(better_neighbors)
                    % insert function values
                    pop.pop_data.neighbor_best_fits(better_neighbors) = ...
                        neighbor_best(better_neighbors);
                    % insert individuals
                    for i = (find(better_neighbors)).'
                        pop.pop_data.neighbor_best_inds(i, :) = ...
                            new_inds(neighbors(i, ind(i)), :);
                    end
                end
                % chop off additional NaN-entry again
                new_fits = new_fits(1:end-1);

                % update the local bests
                new_locals = new_fits < pop.pop_data.local_best_fits;
                pop.pop_data.local_best_fits(new_locals, :) = new_fits(new_locals, :);
                pop.pop_data.local_best_inds(new_locals, :) = new_inds(new_locals, :);

                % update the global best
                if (min(new_fits) < pop.pop_data.global_best_fit)
                    [pop.pop_data.global_best_fit, ind] = min(new_fits);
                    pop.pop_data.global_best_ind = new_inds(ind, :);
                end
                
            % multi objective
            else
                
                % rename for clarity
                new_inds = pop.pop_data.parent_population;
                
                % set of current Pareto solutions
                Paretos = find(pop.pop_data.front_number == 0);

                % rename for clarity
                % NOTE: [NumNeighbors] is equal to the amount of stars or
                % rings,  so using [pop.options.PSO.NumNeighbors] is incorrect!            
                NumNeighbors = numel(pop.pop_data.neighbors(1, :));

                % compute new local and global bests, compute new
                % neighbors and the new best neighbors
                for i = 1:pop.size

                    % new local bests
                    if any(Paretos == i) ||...
                       all(pop.fitnesses(i, :) <= pop.pop_data.local_best_fits(i, :) )
                        % Replace local best when it
                        % - is part of the new Pareto front
                        % - OR it dominates the previous local best
                        pop.pop_data.local_best_fits(i, :) = pop.fitnesses(i, :);
                        pop.pop_data.local_best_inds(i, :) = pop.trans_individuals(i, :);
                    end

                    % New neighbors are [NumNeighbors] particles,
                    % randomly drawn from the Pareto front
                    Paretos_minus_i = Paretos(Paretos ~= i);
                    new_neighs = Paretos_minus_i(round(rand(NumNeighbors,1)*...
                                 (size(Paretos_minus_i,1)-1))+1);
                    pop.pop_data.neighbors(i, :) = new_neighs;

                    % update the neighbor bests
                    [dummy, max_dist] = max(pop.pop_data.crowding_distance(new_neighs));
                    best_neighbor     = new_neighs(max_dist);
                    if all(pop.fitnesses(best_neighbor, :) <= pop.pop_data.neighbor_best_fits(i, :))
                        % replace best neighbor when it
                        % - dominates the previous neighbor best
                        % - is a member of the Pareto front
                        % - is NOT equal to the local best
                        % - has the largest crowding distance of all neighbors
                        pop.pop_data.neighbor_best_fits(i, :) = ...
                            pop.fitnesses(best_neighbor, :);
                        pop.pop_data.neighbor_best_inds(i, :) = ...
                            pop.individuals(best_neighbor, :);
                    end

                end % for
                
                % update the global best
                [maxdist, index] = max(pop.pop_data.crowding_distance(Paretos, :));%#ok
                ffits = pop.fitnesses(Paretos, :);  iinds = pop.trans_individuals(Paretos, :);
                if all(ffits(index(1), :) <= pop.pop_data.global_best_fit)
                    % replace global best when it 
                    % - dominates the previous one
                    % - is a member of the Pareto front, and has
                    % - the largest crowding distance
                    pop.pop_data.global_best_fit = ffits(index(1), :);
                    pop.pop_data.global_best_ind = iinds(index(1), :);
                end % if 
                
            end
            
            % create random matrices
            reparray = ones(1,pop.dimensions);
            r1  = rand(pop.size, 1);  r1 = r1(:, reparray);
            r2  = rand(pop.size, 1);  r2 = r2(:, reparray);
            r3  = rand(pop.size, 1);  r3 = r3(:, reparray);
            
            % update velocities
            pop.pop_data.velocities = ...
                pop.options.PSO.omega  *pop.pop_data.velocities + ...
                pop.options.PSO.eta1   *r1.*(pop.pop_data.neighbor_best_inds - new_inds)+...
                pop.options.PSO.eta2   *r2.*...
                                  bsxfun(@minus, pop.pop_data.global_best_ind, new_inds)+...
                pop.options.PSO.eta3   *r3.*(pop.pop_data.local_best_inds - new_inds);
            
            % limit speeds
            speeds = sqrt(sum(pop.pop_data.velocities.^2, 2));
            velocities_unit_vectors = bsxfun(@rdivide, pop.pop_data.velocities, speeds);
            lower_limit_broken = speeds < pop.pop_data.lower_speed_limit;
            upper_limit_broken = speeds > pop.pop_data.upper_speed_limit;            
            if any(lower_limit_broken)
                pop.pop_data.velocities(lower_limit_broken, :) = ...
                    velocities_unit_vectors(lower_limit_broken, :)*pop.pop_data.lower_speed_limit;
            end
            if any(upper_limit_broken)                
                pop.pop_data.velocities(upper_limit_broken, :) = ...
                    velocities_unit_vectors(upper_limit_broken, :)*pop.pop_data.upper_speed_limit;
            end
            
        end % update algorithm
        
    end % protected/hidden methods
    
end % classdef
