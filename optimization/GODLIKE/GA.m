classdef GA < pop_multi
% GENETIC ALGORITHM
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
        
        % simple constructor: create pop-object
        function pop = GA(varargin)
            % construct pop_single object
            pop = pop@pop_multi(varargin{:});            
        end % constructor
        
        % single iteration
        function iterate(pop)
            
            % the sequence of one iteration depends on whether the
            % optimization is single or multi-objective
            if pop.num_objectives == 1 % single-objective
                pool = pop.tournament_selection_single(pop.size, 2); % binary tournament selection for GA                           
                pop.create_offspring(pool);          % create offspring                             
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
        
        % tournament selection 
        function pool = tournament_selection_single(pop, pool_size, tournament_size)
            
            % initialize mating pool
            pool = zeros(pool_size, 1);            
            % total number of competitors
            rnd_inds = zeros(pool_size*tournament_size,1);
            
            % create random indices outside the loop (faster)
            for i = 1:floor(pool_size*tournament_size/pop.size)
                offset = pop.size*(i-1);
                [dummy, rnds] = sort(rand(pop.size,1));%#ok
                rnd_inds(offset+1:min(end,offset+pop.size), :) = rnds(1:min(end,nnz(~rnd_inds)));
            end
            
            % fill the mating pool
            for i = 1:pool_size                   
                % select [tournament_size] individuals at random
                inds = rnd_inds(1:tournament_size);
                rnd_inds = rnd_inds(tournament_size+1:end);                
                % let them compete according to
                % (xj < yj) if fit(xj) < fit(yj)
                [best, ind] = min(pop.fitnesses(inds));               
                
                % insert the index of the best one in the pool
                pool(i) = inds(ind);              
            end % for
            
        end % function (tournament selection)        
                
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
                                        
            % Generating offspring in GA requires quite many
            % operations. Compared to DE or PSO, it's really
            % quite messy:
            
            % rename some stuff
            Coding     = pop.options.GA.Coding;
            MutProb    = pop.options.GA.MutationProb;
            CrossProb  = pop.options.GA.CrossProb;
            NumBits    = pop.options.GA.NumBits;
            reinit     = pop.options.ReinitRatio;
            if strcmpi(Coding, 'Binary')
                Binary = true;  Real = false;
            else
                Binary = false; Real = true;
            end
            
            % save signs
            signs = sign(parent_pop);
            
            % initialize some arrays that keep track of the signs
            child_signs = zeros(2, pop.dimensions);
            new_signs   = zeros(pop.size, pop.dimensions);
            
            % convert to binary
            if Binary
                % find correct multiplier
                multiplier = 1;
                temp_pop   = round(parent_pop);           % initialize
                parent_pop = abs(parent_pop);             % take absolute value
                temp_pop   = abs(temp_pop);               % take absolute value
                while (max(temp_pop(:)) <= 2^(NumBits)) && ~all(temp_pop(:) == 0)
                    multiplier = multiplier*10;           % adjust multiplier
                    temp_pop   = round(parent_pop*multiplier);% convert to integers
                end
                multiplier = multiplier/10;               % correct multiplier
                temp_pop   = round(parent_pop*multiplier);% convert to integers
                
%                 % see if selected number of bits cannot represent population
%                 if (multiplier == 0.1) && ~all(temp_pop(:) == 0)
%                     error('pop_single:numbits_insufficient', ...
%                         ['Maximum value in population can not be represented by the\n',...
%                         'selected number of bits. Increase ''NumBits'' option, or\n',...
%                         'rescale your problem.']);
%                 end

                % convert each column separately
                bit_representation = false(pool_size, pop.dimensions*NumBits);
                % NOTE: MATLAB's DEC2BIN() should be avoided, as it would
                % be called in a loop and is not a builtin. Moreover, the
                % output of DEC2BIN() is a string, which is pretty
                % inconvenient for the mutation operator. Therefore, do
                % the conversion manually
                for i = 1:pop.dimensions
                    % convert this column to bits
                    bits = temp_pop(:, i);
                    bits = bits*pow2(1-NumBits:0);
                    bits = floor(bits);            % oddly enough, this is much faster
                    bits = bits - fix(bits/2)*2 ;  % than using REM(bits,2)...
                    bits = logical(bits);
                    % append in output matrix
                    bit_representation(:, NumBits*(i-1)+1:NumBits*(i-1)+NumBits) = bits;
                end
                % equate the two
                parent_pop = bit_representation;
                % redefine newpop
                newpop = false(pop.size, pop.dimensions*NumBits);
                % initialize children
                children = false(2, size(parent_pop,2)); 
                % and define convenient conversion array
                % (starts mattering for large population sizes)
                convert_to_dec = repmat(2.^(NumBits-1:-1:0), pop.size, 1);
                
                % convert to array of strings in case of real-representation
            elseif Real
                % do everything in one go with INT2STR()
                % (avoid NUM2STR, as its horrifically slowin a loop)
                % (note that the signs are still included in the array)
                real_representation = int2str(abs(parent_pop)*1e18);
                % convert to array of doubles
                real_representation = real_representation - '0';
                % equate the two
                parent_pop = real_representation;
                % initialize children
                children = zeros(2, size(parent_pop, 2)); 
                % redefine newpop
                newpop = zeros(pop.size, size(parent_pop, 2)); 
            end
            
            % perform crossover
            for i = 1:2:pop.size-1
                
                % select two parents at random
                parent_inds  = round(rand(2,1)*(pool_size-1)+1);
                parents      = parent_pop(parent_inds, :);
                if Binary, parent_signs = signs(parent_inds, :); end
                
                % crossover if a random number is less than [CrossProb]
                % otherwise, just insert the two parents into the new
                % population
                if (rand < CrossProb)
                    
                    % select random crossover point
                    crosspoint = round(rand*(size(parents,2)-1))+1; 
                    
                    % perform crossover
                    children(1, :) = [parents(1,1:crosspoint),parents(2,crosspoint+1:end)];
                    children(2, :) = [parents(2,1:crosspoint),parents(1,crosspoint+1:end)];
                    
                    % also keep track of the signs
                    if Binary
                        index = ceil(crosspoint/NumBits);
                        child_signs(1, :) = [parent_signs(1,1:index),...
                            parent_signs(2,index+1:end)];
                        child_signs(2, :) = [parent_signs(2,1:index),...
                            parent_signs(1,index+1:end)];
                    end
                    
                    % insert children
                    newpop(i:i+1, :)    = children;
                    if Binary, new_signs(i:i+1, :) = child_signs; end
                    
                else
                    newpop(i:i+1, :)    = parents;
                    newfit(i:i+1, :)    = parent_fit(parent_inds, :);
                    if Binary, new_signs(i:i+1, :) = parent_signs; end
                end % if
            end % for
            
            % if the population size is an uneven number, the last entry
            % is still open. Just stick a random parent there
            if mod(pop.size,2)
                index = round(rand*(pool_size-1)+1);
                newpop(end, :)   = parent_pop(index, :);
                newfit(end, :)   = parent_fit(index, :);
                if Binary, new_signs(end,:) = signs(index,:); end
            end
            
            % mutation operator
            mutate = rand(pop.size, size(parent_pop,2)) < MutProb; 
            % If any individual mutates, the function has to be re-evaluated
            newfit(sum(mutate,2)>0,:) = NaN;
            % Binary coding - simply flip bits
            if Binary, newpop(mutate) = ~newpop(mutate); end
            % Real coding - select a new number from [0,9]
            if Real
                % don't mutate spaces
                space_inds = newpop(mutate) == (' '-'0');
                mutate(mutate) = ~space_inds;
                % don't mutate signs
                sign_inds = newpop(mutate) == ('-'-'0');
                mutate(mutate) = ~sign_inds;
                % random new digits
                rnd_inds = round(rand(nnz(mutate),1)*9);
                % convert to strings
                newpop(mutate) = rnd_inds;
            end
            
            % convert back to real numbers
            if Binary
                % initialize
                temp_pop = zeros(pop.size, pop.dimensions);
                % convert columnwise again
                for i = 1:pop.dimensions
                    % convert column to decimal representation
                    temp_pop(:, i) = sum(convert_to_dec.*newpop(:, 1:NumBits), 2);
                    % delete entries
                    newpop = newpop(:, NumBits+1:end);
                end
                % divide by multiplier, and re-assign signs
                temp_pop = temp_pop/multiplier.*new_signs;
                % assign newpop
                newpop = temp_pop;
            elseif Real
                % initialize
                temp_pop = zeros(pop.size, pop.dimensions);
                % assign space character
                space = ' '-'0';
                % append one "space" to the end
                newpop(:, end+1) = space;
                % then convert back to double, column per column
                for i = 1:pop.dimensions
                    % trim leading "spaces"
                    while all(newpop(:,1)==space), newpop = newpop(:, 2:end); end
                    % first find one that does not begin with a "space"
                    non_space = find(newpop(:,1) ~= space, 1);
                    % find indices for the next "space"
                    space_ind = find(newpop(non_space,:) == space, 1);
                    space_ind = space_ind-1;
                    % use power trick forthe conversion
                    powers = 10.^(space_ind-1:-1:0);
                    powers = powers(ones(pop.size,1),:);
                    % remove residual spaces
                    ttemp_pop = newpop(:,1:space_ind);
                    ttemp_pop(ttemp_pop == space) = 0;
                    % insert in final array
                    temp_pop(:, i) = sum(ttemp_pop.*powers,2)/1e18;
                    % adjust newpop
                    newpop = newpop(:, space_ind+1:end);
                end
                % assign newpop
                newpop = signs.*temp_pop;
            end
            
            % reinitialize a small percentage of the population  
            reinit = floor(pop.size * reinit);
            if (reinit > 0)
                % re-initialize
                newpop(1:reinit,:) = pop.lb(1:reinit,:) + ...
                    rand(reinit, pop.dimensions) .* (pop.ub(1:reinit,:)-pop.lb(1:reinit,:));
                % apply sine-transformation
                newpop(1:reinit,:) = real(asin(2*(newpop(1:reinit,:)-pop.lb(1:reinit,:))./...
                    (pop.ub(1:reinit,:)-pop.lb(1:reinit,:))-1));
                % also set corresponding function values to [NaN]
                newfit(1:reinit,:) = NaN;
            end  
            
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
            pop.individuals = pop.wrapperFcn(pop.trans_individuals, 1:pop.size);
            
        end % replace parents
        
    end % protected/hidden methods
    
end % classdef