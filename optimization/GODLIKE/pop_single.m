classdef pop_single < handle
% =insert documentation here=

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

    % public properties
    properties (Access = public)
        algorithm          % type of optimization algorithm used 
        funfcn             % objective function(s)
        confcn             % constraint function(s)
        constrained        % contains type of constraints
        individuals        % the (untransformed) individuals
        fitnesses          % corresponding (penalized) fitnesses        
        size               % population size
        lb                 % lower bounds
        ub                 % upper bounds
        orig_size          % original size of the input
        dimensions         % dimensions
        funevals   = 0;    % number of funct
        iterations = 0;    % iterations so far performed
        options            % options structure (see function [set_options] for info)            
        pop_data           % structure to store intermediate data                           
        eq_indices         % those indices of LU and UB, for which their values are equal
        eq_values          % the corresponding equal values         
        trans_individuals  % transformed individuals 
                           % (sine transformed and without the equal values)
        % contents for single-objective optimization:
        %      pop_data.parent_population
        %      pop_data.offspring_population
        %      pop_data.function_values_parent
        %      pop_data.function_values_offspring        
        %      pop_data.unpenalized_function_values_parent     (for constrained optimization)
        %      pop_data.unpenalized_function_values_offspring  (for constrained optimization)
        %      pop_data.constraint_violations_parent           (for constrained optimization)
        %      pop_data.constraint_violations_offspring        (for constrained optimization)
        
    end
            
    % public methods
    methods (Access = public)
        
        % constructor
        function pop = pop_single(varargin)
                        
            % = = = = = = = = = = = = = = = = = = = = = = = = = = 
            % input is ( new [pop_data] structure, previous [population] object, options )
            % (subsequent call from GODLIKE)
            % = = = = = = = = = = = = = = = = = = = = = = = = = = 
            
            if (nargin == 3) 
                                
                % assign new pop_data structure
                pop.pop_data = varargin{1};
                
                % simply copy previous object
                pop.funfcn      = varargin{2}.funfcn;       pop.iterations = varargin{2}.iterations;
                pop.lb          = varargin{2}.lb;           pop.ub         = varargin{2}.ub;
                pop.funevals    = varargin{2}.funevals;     pop.eq_values  = varargin{2}.eq_values;
                pop.dimensions  = varargin{2}.dimensions;   pop.orig_size  = varargin{2}.orig_size;
                pop.eq_indices  = varargin{2}.eq_indices;   pop.confcn     = varargin{2}.confcn;
                pop.constrained = varargin{2}.constrained;
                                
                % copy individuals and fitnesses
                pop.fitnesses         = pop.pop_data.function_values_parent;                
                pop.trans_individuals = pop.pop_data.parent_population;
                
                % size and options might have changed
                pop.size = size(pop.trans_individuals, 1);%#ok
                pop.options = varargin{3}; 
                
                % set proper algorithm
                pop.algorithm  = pop.options.algorithm;
                
                % replicate [ub], [lb] and the equal values and indices
                % (necessary because of different population size)
                pop.lb = repmat(pop.lb(1, :), pop.size, 1);    
                pop.ub = repmat(pop.ub(1, :), pop.size, 1);
                pop.eq_values  = repmat(pop.eq_values(1, :), pop.size, 1);    
                pop.eq_indices = repmat(pop.eq_indices(1, :), pop.size, 1);
                
                % compute un-transformed individuals
                pop.individuals = pop.wrapperFcn(pop.trans_individuals, 1:pop.size);
                                
                % return
                return
            end
            
            % = = = = = = = = = = = = = = = = = = = = = = = = = = 
            % input is ( funfcn, confcn, popsize, lb, ub, options )
            % (initialization call from GODLIKE)
            % = = = = = = = = = = = = = = = = = = = = = = = = = = 
                         
            % � � � � � � � � � � � � � � � � � � � � � �
            % parse input
            % � � � � � � � � � � � � � � � � � � � � � �
            
            % assign input
            pop.funfcn  = varargin{1};  pop.ub         = varargin{5};
            pop.confcn  = varargin{2};  pop.orig_size  = varargin{6};
            pop.size    = varargin{3};  pop.options    = varargin{7}; 
            pop.lb      = varargin{4};
            
            % set constrained property
            con_in_con = ~all(cellfun(@isempty, pop.confcn));
            con_in_obj = pop.options.ConstraintsInObjectiveFunction > 0;
            if con_in_con && ~con_in_obj
                pop.constrained = 'confcn';
            elseif ~con_in_con && con_in_obj                
                pop.constrained = 'objfun';
            elseif con_in_con && con_in_obj                
                pop.constrained = 'both';
            else
                pop.constrained = [];
            end
            
            % find equal indices and values
            pop.eq_indices = pop.lb == pop.ub;
            pop.eq_values  = pop.lb(pop.eq_indices);
            % the folling command makes this also work for one-dimensional
            % functions, which have no indices equal
            if isempty(pop.eq_values), pop.eq_values = zeros(1,0); end 
                        
            % remove these from [lb] and [ub]
            pop.lb(pop.eq_indices) = [];
            pop.ub(pop.eq_indices) = [];
            
            % replicate [lb], [ub] and the equal indices and values to 
            % facilitate things a bit (and speed it up some more)
            pop.lb = repmat(pop.lb, pop.size, 1);   pop.ub = repmat(pop.ub, pop.size, 1);
            pop.eq_indices = repmat(pop.eq_indices, pop.size, 1); 
            pop.eq_values  = repmat(pop.eq_values, pop.size, 1);
                  
            % NOW set dimensions
            pop.dimensions = size(pop.lb, 2);%#ok
            
            % set optimization algorithm
            pop.algorithm = pop.options.algorithm;
                          
            % � � � � � � � � � � � � � � � � � � � � � � 
            % Initialize population             
            % � � � � � � � � � � � � � � � � � � � � � � 
            
             % initialize population 
            pop.trans_individuals = pop.lb + rand(pop.size, pop.dimensions) .* (pop.ub-pop.lb);
                        
            % apply sine-transformation 
            pop.trans_individuals = real(asin(2*(pop.trans_individuals-pop.lb)./(pop.ub-pop.lb)-1));
            
            % un-transformed individuals 
            pop.individuals = pop.wrapperFcn(pop.trans_individuals, 1:pop.size);

            % insert copy into info structure
            pop.pop_data.parent_population = pop.trans_individuals;            
                        
            % temporarily copy parents to offspring positions 
            pop.pop_data.function_values_offspring = [];
            pop.pop_data.offspring_population      = pop.trans_individuals;
            
            % evaluate function for initial population (parents only)             
            pop.evaluate_function;              
                        
            % copy function values and constraint violations
            pop.fitnesses = pop.pop_data.function_values_offspring; 
            pop.pop_data.function_values_parent = pop.pop_data.function_values_offspring;  
            if ~isempty(pop.constrained)
                pop.pop_data.unpenalized_function_values_parent = pop.pop_data.unpenalized_function_values_offspring;
                pop.pop_data.constraint_violations_parent = pop.pop_data.constraint_violations_offspring;                                
                pop.pop_data.unpenalized_function_values_offspring(:) = inf;
                pop.pop_data.constraint_violations_offspring(:) = inf;
            end
            pop.pop_data.function_values_offspring(:) = inf; 
                                              
        end % function (constructor)  
        
    end % public methods
        
    % protected/hidden methods
    methods (Access = protected, Hidden)
                   
        % wrapper function which includes the equal-valued, and 
        % applies the sine-transformation
        function transformed_individuals = wrapperFcn(pop, input_population, sites)
            % some initializations
            num_sites = nnz(sites);
            transformed_individuals = zeros(num_sites, prod(pop.orig_size));
            % first unapply the sine-transformation
            transformed_individuals(~pop.eq_indices(1:num_sites, :)) = ...
                pop.lb(sites, :) + (pop.ub(sites,:)-pop.lb(sites,:)).*...
                (sin(input_population(sites, :)) + 1)/2;
            % then include the fixed values   
            transformed_individuals(pop.eq_indices(1:num_sites, :)) = pop.eq_values(1:num_sites, :);
        end
        
        % evaluate the objective function(s) correctly
        function evaluate_function(pop)
                   
            % NOTE: suited for both single and multi-objective optimization
            
            % find evaluation sites
            if isempty(pop.pop_data.function_values_offspring)
                sites = (1:pop.size).'; % only in pop-initialization
            else
                sites = ~isfinite(pop.pop_data.function_values_offspring(:, 1));
            end
            
            % number of sites
            num_sites = nnz(sites); iters = 0;
            
            % loop until all function values are finite
            while(any(sites)) && (iters < 25) && (pop.funevals < pop.options.MaxFunEvals)
                                
                % initialize
                fevals = 0; iters = iters + 1;
                
                % include the equal indices and apply sine-transformation
                real_individuals = pop.wrapperFcn(pop.pop_data.offspring_population, sites);
                
                % convert this population to a 3-D cell
                true_pop = reshape(real_individuals.', [pop.orig_size, num_sites]);
                % NOTE: for-loop with JIT is faster than MAT2CELL
                cell_pop = cell(1,1,num_sites);
                for ii = 1:size(true_pop,3), cell_pop{1,1,ii} = true_pop(:,:,ii); end %#ok
                
                % carefully evaluate all functions with cellfun
                arg_out = cell(1, pop.options.ConstraintsInObjectiveFunction+1);
                convals = 0;
                for ii = 1:numel(pop.funfcn)                    
                    try
                        % evaluate function
                        if isempty(pop.constrained) || strcmpi(pop.constrained, 'confcn')  
                            funvals = ...
                                cellfun(@(x)feval(pop.funfcn{ii},x),cell_pop,'uniformoutput',false);
                            % reshape
                            funvals = reshape([funvals{:}],[],num_sites).';
                            % constraint functions might be calculated inside the
                            % objective function(s):
                        elseif any(strcmpi(pop.constrained, {'objfun';'both'}))
                            [arg_out{:}] = ...
                                cellfun(@(x)feval(pop.funfcn{ii},x),cell_pop,'uniformoutput',false);
                            funvals = arg_out{1};
                            new_c_vals   = arg_out{pop.options.ConstraintsInObjectiveFunction};
                            new_c_vals   = new_c_vals(:);
                            new_ceq_vals = arg_out{pop.options.ConstraintsInObjectiveFunction+1};
                            new_ceq_vals = new_ceq_vals(:);
                            % reshape
                            funvals = reshape([funvals{:}],[],num_sites).';
                            % GODLIKE uses the *SUM* of all constraint violations
                            new_c_vals   = cellfun(@(x)sum(x(x>pop.options.TolCon)), new_c_vals);
                            new_ceq_vals = cellfun(@(x)sum(x(abs(x)>pop.options.TolCon)), new_ceq_vals);
                            convals = convals + ...
                                reshape(new_c_vals(:),[],num_sites).' + ...
                                reshape(new_ceq_vals(:),[],num_sites).';
                        end
                        % insert in array
                        pop.pop_data.function_values_offspring(sites, ...
                            ii:ii+size(funvals,2)-1) = funvals; %#ok
                        % keep track of function evaluations
                        fevals = fevals + num_sites;
                        
                        % evaluating objective functions might fail
                    catch userFcn_ME
                        pop_ME = MException('pop_single:function_doesnt_evaluate',...
                            'GODLIKE cannot continue: failure during objective function evaluation.');
                        userFcn_ME = addCause(userFcn_ME, pop_ME);
                        rethrow(userFcn_ME);
                    end
                end
                
                % carefully evaluate all constraint functions with cellfun
                if any(strcmpi(pop.constrained, {'confcn';'both'}))
                    for ii = 1:numel(pop.confcn)
                        try
                            [new_c_vals, new_ceq_vals] = ...
                                cellfun(@(x)feval(pop.confcn{ii},x),cell_pop,'uniformoutput',false);
                            new_c_vals = new_c_vals(:); new_ceq_vals = new_ceq_vals(:);
                            % GODLIKE uses the *SUM* of all constraint violations
                            new_c_vals   = cellfun(@(x)sum(x(x>pop.options.TolCon)), new_c_vals);
                            new_ceq_vals = cellfun(@(x)sum(x(abs(x)>pop.options.TolCon)), new_ceq_vals);
                            convals = convals + ...
                                reshape(new_c_vals(:),[],num_sites).' + ...
                                reshape(new_ceq_vals(:),[],num_sites).';
                            % update fevals counter
                            fevals = fevals + num_sites;
                            % evaluating constraint function might fail
                        catch userFcn_ME
                            pop_ME = MException('pop_single:function_doesnt_evaluate',...
                                'GODLIKE cannot continue: failure during constraint function evaluation.');
                            userFcn_ME = addCause(userFcn_ME, pop_ME);
                            rethrow(userFcn_ME);
                        end
                    end %for                    
                end % if
                
                % penalize and insert constraint violations in pop_data
                if ~isempty(pop.constrained)
                    pop.pop_data.function_values_offspring(sites,:) = pop.penalize(...
                        pop.pop_data.function_values_offspring(sites,:), convals, sites);
                end
                
                % update number of function evaluations
                pop.funevals = pop.funevals + fevals;
                
                % normally, all evaluation sites should now be finite. In
                % some not very well-posed problems, this might however not
                % be the case; some individuals, their function values or 
                % constraint violations might be inf or NaN. In those cases, 
                % we need to re-initialize the corresponding individuals
                
                % Get indices to invalid individuals
                nonfinite_fitnesses    = ~isfinite(pop.pop_data.function_values_offspring);
                nonfinite_individuals  = ~isfinite(pop.pop_data.offspring_population);
                
                % indices to those that need reinitialization
                sites = any(nonfinite_fitnesses,2) | ...
                        any(nonfinite_individuals,2);
                
                % constrained version
                if ~isempty(pop.constrained)
                    nonfinite_constrvalues = ...
                        ~isfinite(pop.pop_data.constraint_violations_offspring);
                    sites = sites | any(nonfinite_constrvalues,2);
                end
                
                % we're done if it's empty
                if ~any(sites), break; end
                         
                % otherwise, re-initialize   
                num_sites = nnz(sites);
                pop.pop_data.offspring_population(sites,:) = pop.lb(sites,:) + ...
                    rand(num_sites, pop.dimensions).*(pop.ub(sites,:)-pop.lb(sites,:));
                % and re-transform
                pop.pop_data.offspring_population(sites,:) = ...
                    real(asin(2*(pop.pop_data.offspring_population(sites,:)-pop.lb(sites,:))...
                    ./(pop.ub(sites,:)-pop.lb(sites,:))-1));
                
            end % evaluation while loop 
            
        end % evaluate function 
        
        % compute penalties and insert constraint violations in pop_data
        % (used for constrained optimizations)
        function funvals = penalize(pop, funvals, convals, sites)
          
            % make sure that the ones that are not violated are zero
            convals(convals <= pop.options.TolCon) = 0;
            
            % insert in pop.pop_data
            pop.pop_data.constraint_violations_offspring(sites, 1) = convals;
            pop.pop_data.unpenalized_function_values_offspring(sites, :) = funvals;
                        
            % assign penalties
            if any(convals(:) > 0)                
                % scaling parameter                
                scale = min(1e16, 1/pop.options.TolCon);                
                % replicate convals
                if (size(funvals,2) > 1), convals = convals(:, ones(size(funvals,2),1)); end%#ok
                
funvals = funvals + scale*convals;
                
%                 % detect which penalties are going to overflow
%                 overflows = convals > 50;
%                 % penalize these with a linear penalty function
%                 funvals(~overflows) = funvals(~overflows) + exp(convals(~overflows)) - 1;                
%                 % Use exponential penalty function for the others
%                 funvals(overflows) = funvals(overflows) +  ...
%                     exp(50) - 1 + scale*convals(overflows);
            end
            
        end % penalize
        
    end % protected/hidden methods
    
end % classdef
