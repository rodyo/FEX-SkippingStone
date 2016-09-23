function [sol, fval, exitflag, output, grad] = ...
        optimize(funfcn, x0, A, b, Aeq, beq, lb, ub, nonlcon, ...
                 strictness, options, algorithm, varargin)
%OPTIMIZE        Optimize general constrained problems using Nelder-Mead algorithm
%
% Usage:
% sol = OPTIMIZE(func, x0) 
% sol = OPTIMIZE(..., x0, lb, ub)
% sol = OPTIMIZE(..., ub, A, b) 
% sol = OPTIMIZE(..., b, Aeq, beq) 
% sol = OPTIMIZE(..., beq, nonlcon) 
% sol = OPTIMIZE(..., nonlcon, strictness) 
% sol = OPTIMIZE(..., strictness, options) 
% sol = OPTIMIZE(..., options, algorithm) 
%
% [sol, fval] = OPTIMIZE(func, ...)
% [sol, fval, exitflag] = OPTIMIZE(func, ...)
% [sol, fval, exitflag, output] = OPTIMIZE(func, ...)
%
% INPUT ARGUMENTS:
%
%  fun, x0, options, varargin - see the help for FMINSEARCH.
%
%  lb - (OPTIONAL) lower bound vector or array, must have the same 
%       size as x0.
%
%       If no lower bounds exist for one of the variables, then
%       supply -inf for that variable.
%
%       If no lower bounds exist at all, then [lb] may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
%  ub - (OPTIONAL) upper bound vector or array, must have the same 
%       size as x0.
%
%       If no upper bounds exist for one of the variables, then
%       supply +inf for that variable.
%
%       If no upper bounds at all, then [ub] may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
%  A, b - (OPTIONAL) Linear inequality constraint array and right
%       hand side vector. (Note: these constraints were chosen to
%       be consistent with those of fmincon.)
%
%       This linear constraint forces the solution vector [x] to 
%       satisfy 
%                               A*x <= b
%
%       Note that in case [x] is a matrix (this is true when [x0] is
%       a matrix), the argument [b] must have corresponding size 
%       [size(A,1) x size(x0,2)], since the same equation is used to 
%       evaluate this constraint. 
%
%  Aeq, beq - (OPTIONAL) Linear equality constraint array and right
%       hand side vector. (Note: these constraints were chosen to
%       be consistent with those of fmincon.)
%
%       This linear constraint forces the solution vector [x] to 
%       satisfy 
%
%                               Aeq*x == beq
%
%       Note that in case [x] is a matrix (this is true when [x0] is
%       a matrix), the argument [beq] must have corresponding size 
%       [size(Aeq,1) x size(x0,2)], since the same equation is used to 
%       evaluate this constraint.
%
%  nonlcon - (OPTIONAL) function handle to general nonlinear constraints,
%       inequality and/or equality constraints.
%
%       [nonlcon] must return two vectors, [c] and [ceq], containing the
%       values for the nonlinear inequality constraints [c] and
%       those for the nonlinear equality constraints [ceq] at [x]. (Note: 
%       these constraints were chosen to be consistent with those of 
%       fmincon.)
%
%       These constraints force the solution to satisfy
%
%               ceq(x)  = 0
%                 c(x) <= 0,
%
%       where [c(x)] and [ceq(x)] are general non-linear functions of [x].
%
% strictness - (OPTIONAL) By default, OPTIMIZE will assume the objective 
%       (and constraint) function(s) can be evaluated at ANY point in 
%       RN-space; the initial estimate does not have to lie in the 
%       feasible region, and intermediate solutions are also allowed to step 
%       outside this area. If your function does not permit such behavior, 
%       set this argument to 'strict'. With 'strict' enabled, the linear 
%       constraints will be satisfied strictly, while the nonlinear 
%       constraints will be satisfied within options.TolCon. 
%
%       If this is also not permissible, use 'superstrict' - then all 
%       nonlinear constraints are also satisfied AT ALL TIMES, and the 
%       objective function is NEVER evaluated outside the feasible area. 
%
%       When using 'strict' or 'superstrict', the initial estimate [x0]
%       MUST be feasible. If it is not feasible, an error is produced
%       before the objective function is ever evaluated. 
%
% algorithm - (OPTIONAL) By default, an embedded version of the 
%       Nelder-Mead algorithm is used. This version is slightly more 
%       robust and internally effecient than the one implemented in 
%       FMINSEARCH. The FMINSEARCH algorithm can still be selected, by
%       setting [algorithm] to 'fminsearch'. 
%
%
% OUTPUT ARGUMENTS:
%
% sol, fval - the solution vector and the corresponding function value,
%       respectively. 
%
% exitflag - (See also the help on FMINSEARCH) A flag that specifies the 
%       reason the algorithm terminated. FMINSEARCH uses only the values
%
%           1    fminsearch converged to a solution x
%           0    Max. # of function evaluations or iterations exceeded
%          -1    Algorithm was terminated by the output function.
%
%       Since OPTIMIZE handles constrained problems, the following 
%       values were added:
%
%           2    All elements in [lb] and [ub] were equal - nothing done
%          -2    Problem is infeasible after the optimization (Some or 
%                any of the constraints are violated at the final 
%                solution).
%          -3    INF or NAN encountered during the optimization. 
%
% output - (See also the help on FMINSEARCH) A structure that contains
%       additional details on the optimization. FMINSEARCH returns
%
%           output.algorithm   Algorithm used           
%           output.iterations  Number of iterations 
%           output.message     Exit message
%
%       in addition, OPTIMIZE returns
%
%         with non-empty [nonlcon] argument: 
%           output.ObjfuncCount     Number of evaluations of the given
%                                   objective function
%           output.ConstrfuncCount  Number of evaluations of the given
%                                   onlinear constraint function
%
%         with no or empty [nonlcon] argument:
%           output.funcCount   Number of function evaluations 
%
%       For constrained problems, the following fields are also present:
%
%           output.constrviolation.lin_ineq
%           output.constrviolation.lin_eq
%           output.constrviolation.nonlin_ineq
%           output.constrviolation.nonlin_eq
%
%       All these fields contain a [M x 2]-cell array. The fist column
%       contains a logical index to the constraints, which is true if the
%       constraint was violated, false if it was satisfied. The second
%       column contains the amount of constraint violation. This amount is
%       equal to zero if the constraint was satisfied within
%       options.TolCon.
%
%
% Notes:
%
%  If options is supplied, then TolX will apply to the transformed
%  variables. All other FMINSEARCH parameters should be unaffected.
%
%  Variables which are constrained by both a lower and an upper
%  bound will use a sin() transformation. Those constrained by
%  only a lower or an upper bound will use a quadratic
%  transformation, and unconstrained variables will be left alone.
%
%  Variables may be fixed by setting their respective bounds equal.
%  In this case, the problem will be reduced in size for FMINSEARCH. 
%
%  If your problem has an EXCLUSIVE (strict) bound constraints which 
%  will not permit evaluation at the bound itself, then you must 
%  provide a slightly offset bound. An example of this is a function 
%  which contains the log of one of its parameters. If you constrain 
%  the variable to have a lower bound of zero, then OPTIMIZE may
%  try to evaluate the function exactly at zero.
%
% EXAMPLES:
%
% rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
% <<Fully unconstrained problem>>
%
% optimize(rosen, [3 3])
% ans =
%    1.0000    1.0000
%
%
% <<lower bound constrained>>
%
% optimize(rosen,[3 3],[2 2],[])
% ans =
%    2.0000    4.0000
%
%
% <<x(2) fixed at 3>>
%
% optimize(rosen,[3 3],[-inf 3],[inf,3])
% ans =
%    1.7314    3.0000
%
%
% <<simple linear inequality: x(1) + x(2) <= 1>>
%
% optimize(rosen,[0 0],[],[],[1 1], 1)
% 
% ans =
%    0.6187    0.3813
% 
%
% <<nonlinear inequality: sqrt(x(1)^2 + x(2)^2) <= 1>>
% <<nonlinear equality  : x(1)^2 + x(2)^3 = 0.5>>
%
% execute this m-file:
%
%   function test_optimize
%        rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
%        options = optimset('TolFun', 1e-8, 'TolX', 1e-8);
%
%        optimize(rosen, [3 3], [],[],[],[],[],[],...
%        @nonlcon, [], options)
%
%   end
%   function [c, ceq] = nonlcon(x)
%        c = norm(x) - 1;
%        ceq = x(1)^2 + x(2)^3 - 0.5;
%   end
%
% ans =
%    0.6513    0.4233
%
%
% Of course, any combination of the above constraints is
% also possible.
%
%
% See also: fminsearch, fminsearchcon, fminsearchbnd, fmincon.


% Current version: v5

% Author info
%{
 Author     : Rody P.S. Oldenhuis
   Delft University of Technology
   E-mail   : oldenhuis@gmail.com

 FMINSEARCHBND, FMINSEARCHCON and most of the
 help for OPTIMIZE witten by
               : John D'Errico
        E-mail : woodchips@rochester.rr.com

OPTIMIZE()
   Last edited : 12/Oct/2009
 Last Uploaded : 05/Aug/2009

 [ Please report bugs to oldenhuis@gmail.com ]
%}

% Known problems
%{
    - 'iter-detailed' doesn't work in NelderMead.
      It is useful to show constraint violation etc. when
      iter-detailed is selected.
    - randomized initial population might not be the best choice
      for the global routine...
    - grad_Z() transofrmation function doesn't work correctly with 
      one or more fixed values
%}

%CHANGELOG
%{

Nov 2, 2010
 v5:  - Made output functions work correctly in the 
        global routine. 
      - fixed a bug: fixed variables were inserted 
        incorrectly in the transformation function, when 
        using NelderMead(), resulting in an error.
      - If at any one iteration in NelderMead al function 
        values and/or x-values are non-finite, the algorithm
        terminates. This is different from how fminsearch()
        handles such a case; fminsearch keeps on shrinking 
        the simplex, which is pretty useless most of
        the time.
      - made the overall structure more readible
      - Added option 'popsize' for the global routine
      - Fixed problem with outputFcn, in that 
        optimvalues.fval showed the
        PENALIZED function value, which should of 
        course be the UNpenalized one
      - Fixed problem with globalized version destroying 
        the ConstraintsInObjectiveFunction field in the 
        options structure (apparently, options = 
        optimset(options, 'parameter', value) removes 
        all non-standard fields...)
 
Oct 12, 2009
 v4: - Now, OPTIMIZE() also supports FMINLBFGS(), 
       a limited-memory, BFGS quasi-newton 
       optimizer, also available on the file 
       exchange. This is a major change, because
       it requires more options and, more 
       importantly, derivatives for the penalized
       constraint functions.
     - Penalties are now scaled with the objective 
       function value
     - found lots of small bugs:
        · text formatting in output.message was all 
          messed up for constrained problems
        · ConstrFuncCount wasn't updated; it always
          resulted in OPTIMIZE() reporting only 2  
          evaluations of the constraint function
        · constrained global searches had a lot of 
          problems with finalizing their solutions; 
          these were mostly related to the proper
          formatting of [output].
     - fixed another small bug; [lb] and [ub] were 
       swapped during initialization, resulting in 
       an error when [lb] was given but [ub] was 
       omitted. 
     - Changed the default algorithm back to 
       FMINSEARCH. This is of course the least buggy,
       and most well-known to everyone. The internal 
       NelderMead has now been published separately
       (see MATLAB file exchange; NelderMead)
     - Corrected small mistake in the penalty 
       function for linear equality constraints; 
       an ABS() was omitted, so that the sum of the 
       constraint violation could result in lower 
       penalties than actually deserved.
     

 Aug 6, 2009
 v3: - removed dependency on the OPTIMSET() from the 
       optimization toolbox (TolCon is not part of 
       the diet-OPTIMSET most people have)
     - Included basic global optimization routine. If
       [x0] is omitted, but [ub] and [lb] are given, 
       a number of initial values will be optimized, 
       which are randomly generated in the boundaries 
       [lb] and [ub].
     - added one exitflag (-3). This is the exitflag
       for when the globalized algorithm has found 
       nothing else than INF or NAN function values.


 Aug 5, 2009
 v2: - Included NelderMead algorithm, a slightly more
       robust and internally efficient method than 
       FMINSEARCH. Of course, FMINSEARCH can still be
       selected through yet another argument...
     - one-dimensional functions could not be 
       optimized because of a bad reference to 
       elements of [lb] and [ub]; problem fixed. 
     - I forgot that [strictness] and [options] can 
       also be empty; included an ISEMPTY() check. 
     - fixed a problem with the outputFcn; this 
       function was still evaluated under the 
       assumption that [x0] is always a vector. 
       Inserted a reshaping operation. 


 Aug 1, 2009
 v1: - x0 can now be a matrix, just as in FMINSEARCH
     - minor bug fixed with the strictness setting: 
       with multiple nonlinear constraints, an ANY()
       was necessary, which was not used. 
     
  
 May 24, 2009.
 v0: - first release; simple & elegant.

%}
    
    %% initialize
    
    % process input
    narg = nargin;
    error(nargchk(2, inf, narg));
    if (narg < 12) || isempty(algorithm ), algorithm = 'fminsearch'; end   
    if (narg < 11) || isempty(options   ), options = optimset; end   
    if (narg < 10) || isempty(strictness), strictness = 'loose'; end   
    if (narg <  9), nonlcon = []; end  
    if (narg <  8), ub  = []; end
    if (narg <  7), lb  = []; end 
    if (narg <  6), beq = []; end
    if (narg <  5), Aeq = []; end
    if (narg <  4), b   = []; end
    if (narg <  3), A   = []; end
          
    % Are the constraints computed in the objective function, or are
    % they separate functions?
    if isfield(options, 'ConstraintsInObjectiveFunction')
         confcn_in_objfcn = options.ConstraintsInObjectiveFunction;
    else confcn_in_objfcn = false;
    end
                 
    % get tolerance on constraints
    % NOTE: don't use optimget- it destroys non-standard fields
    if isfield(options, 'TolCon') && ~isempty(options.TolCon)
         tolCon = options.TolCon; % user-provided value
    else tolCon = 1e-6;           % default value    
    end
    
    % set some logicals for easier reading
    have_nonlconFcn    = ~isempty(nonlcon);
    have_linineqconFcn = ~isempty(A) && ~isempty(b);
    have_lineqconFcn   = ~isempty(Aeq) && ~isempty(beq);
    % see if we have to keep track of algorithm data
    create_output = (nargout == 4);
    % check if gradients are required
    compute_grad = strcmpi(algorithm,'fminlbfgs') && ...
        isfield(options,'GradObj') && strcmpi(options.GradObj,'on');
    compute_grad_con = strcmpi(algorithm,'fminlbfgs') && ...
        isfield(options,'GradConstr') && strcmpi(options.GradConstr,'on');
    % do we display anything? 
    do_display = isfield(options, 'Display') && ~isempty(options.Display) && ...
        ~strcmpi(options.Display,'off');
    do_extended_display = do_display && strcmpi(options.Display,'iter-detailed');
    % global optimization or local? 
    do_global_opt = isempty(x0);
    % do we have an output function? 
    have_outputFcn = isfield(options, 'OutputFcn') && ~isempty(options.OutputFcn);

    % no x0 given means optimize globally.
    if do_global_opt
        [sol, fval, exitflag, output] = optimize_globally(varargin{:}); return; end
            
    % make copy of UNpenalized function value
    UPfval = inf;    
    
    % initialize & define constants    
    superstrict = false;      % initially, don't use superstrict setting    
    exp50   = exp(50);        % maximum penalty
    N0      = numel(x0);      % variable to check sizes etc.
    Nzero   = zeros(N0, 1);   % often-used zero-matrix    
    grad    = [];             % initially, nothing for gradient
    sumAT   = repmat(sum(A,1).',1,size(b,2));     % column-sum of [A]  , transposed and replicated
    sumAeqT = repmat(sum(Aeq,1).',1,size(beq,2)); % column-sum of [Aeq], transposed and replicated
    
    % initialize output
    output = [];
    if create_output
        % fields always present
        output.iterations = 0;
        output.algorithm  = '';
        output.message    = 'Initializing optimization...';
        % fields present depending on the presence of nonlcon
        if ~have_nonlconFcn
            output.funcCount = 0;
        else
            output.ObjfuncCount = 0;
            output.ConstrfuncCount = 1; % one evaluation in check_input()
            output.constrviolation.nonlin_eq   = cell(2,1);
            output.constrviolation.nonlin_ineq = cell(2,1);
        end
        % fields present depending on the presence of A or Aeq
        if have_linineqconFcn
            output.constrviolation.lin_ineq = cell(2,1); end
        if have_lineqconFcn
            output.constrviolation.lin_eq = cell(2,1); end
    end
    
    % check for an output function. If there is any, 
    % use wrapper function to call it with un-transformed variable
    if have_outputFcn
        OutputFcn = options.OutputFcn;
        options.OutputFcn = @OutputFcn_wrapper;
    end    
    
    % adjust bounds when they are empty
    if isempty(lb), lb = -inf(size(x0)); end
    if isempty(ub), ub = +inf(size(x0)); end
       
    % check the user-provided input with nested function check_input
    check_input;
        
    % force everything to be column vector
    new_x = x0;      ub = ub(:);     
    x0 = x0(:);      lb = lb(:);     x0one = ones(size(x0));           
    
    % replicate lb or ub when they are scalars, and x0 is not
    if isscalar(lb) && (N0 ~= 1), lb = lb*x0one; end
    if isscalar(ub) && (N0 ~= 1), ub = ub*x0one; end    
    
    % determine the type of bounds
    nf_lb   = ~isfinite(lb);                nf_ub   = ~isfinite(ub); 
    fix_var =  lb == ub;                    lb_only = ~nf_lb &  nf_ub & ~fix_var;
    ub_only =  nf_lb & ~nf_ub & ~fix_var;   unconst =  nf_lb &  nf_ub & ~fix_var;
    lb_ub   = ~nf_lb & ~nf_ub & ~fix_var;    
    
    % if all variables are fixed, simply return
    if length(lb(fix_var)) == N0
        sol  = reshape(lb,size(new_x)); 
        fval = funfcn(sol); 
        exitflag = 2; 
        if create_output
            output.iterations = 0;
            output.message = 'Lower and upper bound were set equal - nothing to do. ';
            if ~have_nonlconFcn
                 output.funcCount = 1;
            else output.ObjfuncCount = 1;
            end
        end
        [output, exitflag] = finalize(lb, output, exitflag);
        if create_output && exitflag ~= -2
            output.message = sprintf(...
                '%s\nFortunately, the solution is feasible using OPTIONS.TolCon of %1.6f.',...
                output.message, tolCon);
        end
        return
    end    
    
    %% procedure
    
    % force the initial estimate inside the given bounds
    x0(x0 < lb) = lb(x0 < lb);  x0(x0 > ub) = ub(x0 > ub);
        
    % transform initial estimate to its unconstrained counterpart
    xin = x0;                                        % fixed and unconstrained variables   
    xin(lb_only) = sqrt(x0(lb_only) - lb(lb_only));  % lower bounds only   
    xin(ub_only) = sqrt(ub(ub_only) - x0(ub_only));  % upper bounds only
    xin(lb_ub)   = real(asin( 2*(x0(lb_ub) - lb(lb_ub))./ ...
                     (ub(lb_ub) - lb(lb_ub)) - 1));  % both upper and lower bounds
    xin(fix_var) = [];
    
    % some more often-used matrices
    None   = ones(numel(xin)+1,1);  
    Np1zero = zeros(N0, numel(xin)+1);
    
    % optimize the problem
    switch lower(algorithm)
        
        % MATLAB's own derivative-free Nelder-Mead algorithm (FMINSEARCH())
        case 'fminsearch'
            % optimize
            [presol, fval, exitflag, output_a] = ...
                fminsearch(@funfcnT, xin, options, varargin{:});                          
            % transform solution back to original (bounded) variables
            sol = new_x;    sol(:) = X(presol); % with the same size as the original x0
            % and evaluate function once more to get unconstrained values
            fval(:) = funfcn(sol);
          
        % a slightly improved version of the Nelder-Mead algorithm
        % (basically the same as FMINSEARCH, but then a bit more internally efficient)
        case 'neldermead'
            % optimize
            [presol, fval, exitflag, output_a] = ...
                NelderMead(@funfcnT, xin, options, varargin{:});
            % transform solution back to original (bounded) variables
            sol = new_x;    sol(:) = X(presol); % with the same size as the original x0
            % and evaluate function once more to get unconstrained values
            fval(:) = funfcn(sol);
                        
        % Steepest descent or Quasi-Newton (limited-memory) BFGS 
        % (both using gradient) FMINLBFGS(), by Dirk-Jan Kroon 
        case 'fminlbfgs'
            % optimize
            [presol, fval, exitflag, output_a, grad] = ...
                fminlbfgs(@funfcnT, xin, options, varargin{:});
            % transform solution back to original (bounded) variables
            sol = new_x;    sol(:) = X(presol); % with the same size as the original x0  
            % and evaluate function some more to get unconstrained values
            if compute_grad
                % function value and gradient
                [fval, grad] = funfcn(sol);                
            else
                % only function value
                fval(:) = funfcn(sol);    grad = zeros(size(grad));
                % compute gradient with simple central differences                
                % NOTE: re-using [sol] in stead of copying to
                % a new variable is more memory efficient
                perturb = 1e-6;
                for k = 1:numel(sol)
                    % add & evaluate
                    sol(k) = sol(k) + perturb;     fval_plus = funfcn(sol);
                    % subtract & evaluate
                    sol(k) = sol(k) - 2*perturb;   fval_minus = funfcn(sol);
                    % compute central differences
                    grad(k) = (fval_plus - fval_minus)/2/perturb;
                    % correct [sol] 
                    sol = sol + perturb;   
                    % add one extra function evaluation (another one is
                    % added below)
                    if create_output
                        output_a.funcCount = output_a.funcCount + 1; end
                end
            end
            
    end % switch (algorithm)
    
    % copy appropriate fields to the output structure
    if create_output
        output.message    = output_a.message;
        output.algorithm  = output_a.algorithm;
        output.iterations = output_a.iterations;        
        if ~have_nonlconFcn
             output.funcCount = output_a.funcCount + 1;
        else output.ObjfuncCount = output_a.funcCount + 1; 
        end
    end

    % append constraint violations to the output structure, and change the
    % exitflag accordingly
    [output, exitflag] = finalize(sol, output, exitflag);
        
    %% NESTED FUNCTIONS (THE ACTUAL WORK)
    
    % check user provided input
    function check_input
        
        % dimensions & weird input
        if (numel(lb) ~= N0 && ~isscalar(lb)) || (numel(ub) ~= N0 && ~isscalar(ub))
            error('optimize:lb_ub_incompatible_size',...
                'Size of either [lb] or [ub] incompatible with size of [x0].')
        end
        if ~isempty(A) && isempty(b)
            warning('optimize:Aeq_but_not_beq', ...
                ['I received the matrix [A], but you omitted the corresponding vector [b].',...
                '\nI`ll assume a zero-vector for [b]...'])
            b = zeros(size(A,1), size(x0,2));
        end
        if ~isempty(Aeq) && isempty(beq)
            warning('optimize:Aeq_but_not_beq', ...
                ['I received the matrix [Aeq], but you omitted the corresponding vector [beq].',...
                '\nI`ll assume a zero-vector for [beq]...'])
            beq = zeros(size(Aeq,1), size(x0,2));
        end
        if isempty(Aeq) && ~isempty(beq)
            warning('optimize:beq_but_not_Aeq', ...
                ['I received the vector [beq], but you omitted the corresponding matrix [Aeq].',...
                '\nI`ll ignore the given [beq]...'])
            beq = [];
        end
        if isempty(A) && ~isempty(b)
            warning('optimize:b_but_not_A', ...
                ['I received the vector [b], but you omitted the corresponding matrix [A].',...
                '\nI`ll ignore the given [b]...'])
            b = [];
        end
        if have_linineqconFcn && size(b,1)~=size(A,1)
            error('optimize:b_incompatible_with_A',...
                'The size of [b] is incompatible with that of [A].')
        end
        if have_lineqconFcn && size(beq,1)~=size(Aeq,1)
            error('optimize:b_incompatible_with_A',...
                'The size of [beq] is incompatible with that of [Aeq].')
        end
        if ~isvector(x0) && ~isempty(A) && (size(A,2) ~= size(x0,1))
            error('optimize:A_incompatible_size',...
                'Linear constraint matrix [A] has incompatible size for given [x0].')
        end
        if ~isvector(x0) && ~isempty(Aeq) && (size(Aeq,2) ~= size(x0,1))
            error('optimize:Aeq_incompatible_size',...
                'Linear constraint matrix [Aeq] has incompatible size for given [x0].')
        end
        if ~isempty(b) && size(b,2)~=size(x0,2)
            error('optimize:x0_vector_but_not_b',...
                  'Given linear constraint vector [b] has incompatible size with given [x0].')
        end
        if ~isempty(beq) && size(beq,2)~=size(x0,2)
            error('optimize:x0_vector_but_not_beq',...
                  'Given linear constraint vector [beq] has incompatible size with given [x0].')
        end                
          
        % functions are not function handles
        if ~isa(funfcn, 'function_handle')
            error('optimize:func_not_a_function',...
                  'Objective function must be given as a function handle.')
        end
        if ~isempty(nonlcon) && ~isa(nonlcon, 'function_handle')
            error('optimize:nonlcon_not_a_function', ...
                 'non-linear constraint function must be a function handle.')
        end
        
        % check the given algorithm
        if ~isempty(algorithm) && ...
           ~any(strcmpi(algorithm, {'fminlbfgs'; 'fminsearch'; 'NelderMead'}))
            error('optimize:invalid_algorithm',...
                ['Invalid algorithm selected. ',...
                'Valid options are ''FMINSEARCH'', ''NELDERMEAD'' or ''FMINLBFGS''.']);
        end        
        if ~isempty(algorithm) && strcmpi(algorithm,'fminlbfgs') && isempty(which('fminlbfgs'))
            error('optimize:fminlbfgs_not_present',...
                'The function FMINLBFGS is not present in the current MATLAB path.');
        end
        
        % evaluate the non-linear constraint function on 
        % the initial value, to perform initial checks
        grad_c = [];   grad_ceq = [];
        if have_nonlconFcn
            if isfield(options, 'GradConstr') && strcmpi(options.GradConstr, 'on')
                [c, ceq, grad_ceq, grad_ceq] = nonlcon(x0);%#ok
            else
                [c, ceq] = nonlcon(x0);%#ok
            end
        end
        
        % check sizes of derivatives
        if ~isempty(grad_c) && (size(grad_c,2) ~= numel(x0)) && (size(grad_c,1) ~= numel(x0))
            error('optimize:grad_c_incorrect_size',...
                ['The matrix of gradients of the non-linear in-equality constraints\n',...
                 'must have one of its dimensions equal to the number of elements in [x].']);            
        end
        if ~isempty(grad_ceq) && (size(grad_ceq,2) ~= numel(x0)) && (size(grad_ceq,1) ~= numel(x0))
             error('optimize:grad_ceq_incorrect_size',...
                ['The matrix of gradients of the non-linear equality constraints\n',...
                 'must have one of its dimension equal to  the number of elements in [x].']);
        end        
        
        % test the feasibility of the initial solution (when strict or
        % superstrict behavior has been enabled)        
        if strcmpi(strictness, 'strict') || strcmpi(strictness, 'superstrict')
            superstrict = strcmpi(strictness, 'superstrict');           
            if ~isempty(A) && any(any(A*x0 > b))
                error('optimize:x0_doesnt_satisfy_linear_ineq', ...
                    ['Initial estimate does not satisfy linear inequality.', ...
                    '\nPlease provide an initial estimate inside the feasible region.']);
            end
            if ~isempty(Aeq) && any(any(Aeq*x0 ~= beq))
                error('optimize:x0_doesnt_satisfy_linear_eq', ...
                    ['Initial estimate does not satisfy linear equality.', ...
                    '\nPlease provide an initial estimate inside the feasible region.']);
            end
            if have_nonlconFcn
                % check [c]
                if ~isempty(c) && any(c(:) > ~superstrict*tolCon) 
                    error('optimize:x0_doesnt_satisfy_nonlinear_ineq', ...
                        ['Initial estimate does not satisfy nonlinear inequality.', ...
                        '\nPlease provide an initial estimate inside the feasible region.']);
                end
                % check [ceq]
                if ~isempty(ceq) && any(abs(ceq(:)) >= ~superstrict*tolCon) 
                    error('optimize:x0_doesnt_satisfy_nonlinear_eq', ...
                        ['Initial estimate does not satisfy nonlinear equality.', ...
                        '\nPlease provide an initial estimate inside the feasible region.']);
                end
            end
        end    
        
    end % check_input
    
    % create transformed variable X to conform to upper and lower bounds
    function Z = X(x)        
        % initialize 
        if (size(x,2) == 1)
            Z = Nzero;    rep = 1;
        else
            Z = Np1zero;  rep = None;
        end        
        % first insert fixed values
        y = x0one(:, rep);
        y( fix_var,:) = lb(fix_var,rep);
        y(~fix_var,:) = x;
        x = y;         
        % and transform
        Z(lb_only, :) = lb(lb_only, rep) + x(lb_only, :).^2;
        Z(ub_only, :) = ub(ub_only, rep) - x(ub_only, :).^2;
        Z(fix_var, :) = lb(fix_var, rep);
        Z(unconst, :) = x(unconst, :);
        Z(lb_ub, :)   = lb(lb_ub, rep) + (ub(lb_ub, rep)-lb(lb_ub, rep)) .* ...
            (sin(x(lb_ub, :)) + 1)/2;        
    end % X
    
    % derivatives of transformed X 
    function grad_Z = gradX(x, grad_x)
        grad_Z             = grad_x;
        grad_Z(lb_only, :) = +2*grad_x(lb_only, :).*x(lb_only, :);
        grad_Z(ub_only, :) = -2*grad_x(ub_only, :).*x(ub_only, :);
        grad_Z(unconst, :) = grad_x(unconst, :);
        grad_Z(lb_ub, :)   = grad_x(lb_ub,:).*(ub(lb_ub,:)-lb(lb_ub,:)).*cos(x(lb_ub,:))/2;
    end % grad_Z
    
    % create penalized function. Penalize with linear penalty function if 
    % violation is severe, otherwise, use exponential penalty. If the
    % 'strict' option has been set, check the constraints, and return INF
    % if any of them are violated.
    function [P_fval, grad_val] = funfcnP(x, varargin)
        
        % initialize function value
        if (size(x,2) == 1), P_fval = 0; else P_fval = None.'-1; end
                        
        % initialize x_new array
        x_new = new_x;
        
        % initialize gradient when needed
        if compute_grad, grad_val = zeros(size(x)); end
        
        % evaluate every column in x
        for i = 1:size(x, 2)
            
            % reshape x, so it has the same size as the given x0
            x_new(:) = x(:, i);
            
            % initialize
            c = [];  ceq = []; grad_c = []; grad_ceq = [];
            
            % evaluate the objective function, taking care that also 
            % a gradient may be supplied
            if compute_grad                
                if ~confcn_in_objfcn
                    [obj_fval, obj_gradient] = funfcn(x_new, varargin{:});                    
                else
                    if compute_grad_con
                         arg_out = cell(1, confcn_in_objfcn+1);
                    else arg_out = cell(1, confcn_in_objfcn+3);
                    end
                    [arg_out{:}] = funfcn(x_new, varargin{:});                    
                    obj_fval     = arg_out{1};
                    obj_gradient = arg_out{2};
                    c            = arg_out{confcn_in_objfcn+0};
                    ceq          = arg_out{confcn_in_objfcn+1};
                    if compute_grad_con
                        grad_c   = arg_out{confcn_in_objfcn+2};
                        grad_ceq = arg_out{confcn_in_objfcn+3};
                    else
                        grad_c   = ''; % use strings to distinguish them later on
                        grad_ceq = '';
                    end
                end                
            else
                if ~confcn_in_objfcn
                    obj_fval = funfcn(x_new, varargin{:});     
                else                    
                    arg_out = cell(1, confcn_in_objfcn+1);
                    [arg_out{:}] = funfcn(x_new, varargin{:});
                    obj_fval     = arg_out{1};
                    c            = arg_out{confcn_in_objfcn+0};
                    ceq          = arg_out{confcn_in_objfcn+1};
                    grad_c       = ''; 
                    grad_ceq     = ''; % use strings to distinguish them later on
                end
                obj_gradient = 1; % needed in penalty function                
            end  
            % amount of function evaluations is already kept track of by
            % FMINSEARCH() or FMINLBFGS()
            
            % make global copy
            UPfval = obj_fval; 
            
            % initially, we are optimistic
            linear_eq_Penalty   = 0;    linear_ineq_Penalty_grad = 0;        
            linear_ineq_Penalty = 0;    linear_eq_Penalty_grad   = 0;
            nonlin_eq_Penalty   = 0;    nonlin_eq_Penalty_grad   = 0;            
            nonlin_ineq_Penalty = 0;    nonlin_ineq_Penalty_grad = 0;  
            
            % Penalize the linear equality constraint violation 
            % required: Aeq*x = beq   
            if have_lineqconFcn
                lin_eq = Aeq*x_new - beq;     
                sumlin_eq = sum(abs(lin_eq(:)));
                if strcmpi(strictness, 'strict') && any(abs(lin_eq) > 0)
                    P_fval = inf; grad_val = inf; return; end
                % compute penalties
                linear_eq_Penalty = Penalize(sumlin_eq, []);
                % also compute derivatives 
                % (NOTE: since the sum of the ABSOLUTE values is used
                % here, the signs are important!)
                if compute_grad
                    linear_eq_Penalty_grad = ...
                        Penalize(sign(lin_eq).*sumAeqT, sumlin_eq); 
                end
            end
                                                        
            % Penalize the linear inequality constraint violation 
            % required: A*x <= b
            if have_linineqconFcn
                lin_ineq = A*x_new - b;                     
                lin_ineq(lin_ineq <= 0) = 0;
                sumlin_ineq = sum(lin_ineq(:)); 
                if strcmpi(strictness, 'strict') && any(lin_ineq > 0)
                    P_fval = inf; grad_val = inf; return; end
                % compute penalties
                linear_ineq_Penalty = Penalize(sumlin_ineq, []);
                % also compute derivatives
                if compute_grad, linear_ineq_Penalty_grad = ...
                        Penalize(sumAT, sumlin_ineq); end                
            end
            
            % Penalize the non-linear constraint violations
            % required: ceq = 0 and c <= 0
            if have_nonlconFcn && ~confcn_in_objfcn
                % initialize as characters, to distinguish them later on; 
                % derivatives may be empty, inf, or NaN as returned from [nonlcon]
                grad_c = ''; grad_ceq = ''; 
                % gradients are given explicitly by [nonlcon]
                if compute_grad_con
                    [c, ceq, grad_c, grad_ceq] = feval(nonlcon, x_new);
                % the gradients are not given by [nonlcon]; they have to 
                % be computed by central differences
                else                    
                    [c, ceq] = feval(nonlcon, x_new);
                    % central-difference derivatives are computed later; 
                    % the strictness setting might make computing it here
                    % unneccecary
                end
                % keep track of number of evaluations made
                if create_output
                    output.ConstrfuncCount = output.ConstrfuncCount + 1; end                
            end
            
            % force grad_c] and [grad_ceq] to be of proper size
            if ~isempty(grad_c)
                grad_c = reshape(grad_c(:), numel(c), numel(x_new)); end
            if ~isempty(grad_ceq)
                grad_ceq = reshape(grad_ceq(:), numel(ceq), numel(x_new)); end

            % process non-linear inequality constraints
            if ~isempty(c)
                % Force it to be a column vector
                c = c(:); 
                % check for strictness setting
                if (strcmpi(strictness, 'strict') || ...
                        strcmpi(strictness, 'superstrict')) &&...
                        any(c > ~superstrict*tolCon)
                    P_fval = inf; grad_val = inf; return
                end  
                % sum the violated constraints
                violated_c = c > tolCon;
                sumc = sum(c(violated_c));
                % compute penalty
                nonlin_ineq_Penalty = Penalize(sumc, []);
            end

            % process non-linear equality constraints
            if ~isempty(ceq)   
                % use the absolute values, but save the signs for the
                % derivatives
                signceq = repmat(sign(ceq), 1, numel(x_new)); 
                ceq = abs(ceq(:)); 
                % check for strictness setting
                if (strcmpi(strictness, 'strict') || ...
                        strcmpi(strictness, 'superstrict')) &&...
                        any(ceq >= ~superstrict*tolCon)
                    P_fval = inf; grad_val = inf; return
                end  
                % sum the violated constraints
                violated_ceq = (ceq >= tolCon); 
                sumceq = sum(ceq(violated_ceq));
                % compute penalty 
                nonlin_eq_Penalty = Penalize(sumceq, []);                    
            end

            % compute derivatives with central-differences of non-linear constraints
            if compute_grad && ischar(grad_c) && ischar(grad_ceq)
                % re-initialize [grad_c] and [grad_ceq]; they're empty (char)'s now 
                % to force OPTIMIZE() to enter this IF-block
                grad_c   = zeros(numel(c), numel(x_new)); 
                grad_ceq = zeros(numel(ceq), numel(x_new));
                % compute central differences
                % NOTE: re-using [x_new] in stead of copying to 
                % a new variable is more memory efficient
                perturb = 1e-6; % default finite difference
                for j = 1:numel(x_new)                        
                    % add perturbation & evaluate
                    x_new(j) = x_new(j) + perturb;   
                    [c_plus , ceq_plus ] = feval(nonlcon, x_new);
                    % subtract perturbation & evaluate
                    x_new(j) = x_new(j) - 2*perturb; 
                    [c_minus, ceq_minus] = feval(nonlcon, x_new);
                    % reset [x_new]
                    x_new(j) = x_new(j) + perturb;                        
                    % compute central differences
                    if ~isempty(c_plus)
                        grad_c(:, j) = (c_plus - c_minus)/2/perturb; end
                    if ~isempty(ceq_plus)
                        grad_ceq(:, j) = (ceq_plus - ceq_minus)/2/perturb; end
                    % keep track of number of evaluations made
                    if create_output
                        output.ConstrfuncCount = output.ConstrfuncCount + 2; end
                end
            end  
            
            % add derivatives of non-linear equality constraint function
            if compute_grad && ~isempty(c)
                % first, remove those that satisfy the constraints
                grad_c = grad_c(violated_c, :);  
                % compute derivatives of penalty functions
                nonlin_ineq_Penalty_grad = ...
                    Penalize(sum(grad_c,1), nonlin_ineq_Penalty);                                        
            end

            % add derivatives of non-linear equality constraint function
            if compute_grad && ~isempty(ceq)
                % first, remove those that satisfy the constraints
                grad_ceq = grad_ceq(violated_ceq, :);
                % compute derivatives of penalty functions
                % (NOTE: since the sum of the ABSOLUTE values is used
                % here, the signs are important!)
                nonlin_eq_Penalty_grad = ...
                    Penalize(sum(signceq.*grad_ceq,1), nonlin_eq_Penalty); 
            end                
            
            % return penalized function value
            P_fval(i) = obj_fval + linear_eq_Penalty + linear_ineq_Penalty + ...
                nonlin_eq_Penalty + nonlin_ineq_Penalty; %#ok MLINT is wrong here...
            
             % return penalized derivatives of constraints             
            grad_val(:, i) = obj_gradient(:) + linear_eq_Penalty_grad(:) + ...
              linear_ineq_Penalty_grad(:) + nonlin_eq_Penalty_grad(:) + nonlin_ineq_Penalty_grad(:);
            
        end
        
        % compute deserved penalties, and derivatives thereof
        % (doubly-nested function)
        function varargout = Penalize(violation, prev_penalty) 
            % scaling parameter
            if isfield(options, 'TolCon') && ~isempty(options.TolCon)
                scale = min(1e16, 1/options.TolCon);
            else
                scale = 1;
            end
            % Compute penalty function
            if isempty(prev_penalty)
                % linear penalty to avoid overflow
                if scale*violation > 50
                    penalty = exp50*(1 + scale*violation)-1;
                % exponential penalty otherwise
                else
                    penalty = exp(scale*violation) - 1;
                end
                varargout{1} = penalty;
            % Compute derivative of penalty function
            else
                % derivative of linear penalty function
                if prev_penalty > exp50
                    derivative = scale*violation*exp50;
                % derivative of exponential penalty function
                else
                    derivative = scale*(prev_penalty+1).*violation;
                end
                varargout{1} = derivative;
            end
        end % Penalize
        
    end % funfcnP
    
    % define the transformed & penalized function    
    function varargout = funfcnT(x, varargin)
        % compute transformed variable
        XT = X(x);
        % with gradient
        if compute_grad
            [varargout{1}, grad_val] = funfcnP(XT, varargin{:});            
            % transform gradient and output
            varargout{2} = gradX(x, grad_val);
        % without gradient
        else
            varargout{1} = funfcnP(XT, varargin{:});
        end 
    end % funfcnT
    
    % simple wrapper function for output functions; these need to be 
    % evaluated with the untransformed variables
    function stop = OutputFcn_wrapper(x, optimvalues, state)
        % transform x        
        x_new = new_x;  x_new(:) = X(x);
        % unpenalized function value
        optimvalues.fval = UPfval;
        % evaluate output function
        stop = OutputFcn(x_new, optimvalues, state);
    end % OutputFcn_wrapper
    
    % finalize the output
    function [output, exitflag] = finalize(x, output, exitflag)
        
        % reshape x so it has the same size as x0
        x_new = new_x; x_new(:) = x;
            
        % compute violations (needed in both display and output structure)
                  
        % initialiy we're optimistic
        is_violated   = false;
        max_violation = 0;
        
        % add proper [constrviolation] field
        if have_linineqconFcn
            Ax        = A*x_new;
            violated  = Ax >= b + tolCon;
            violation = Ax - b;
            violation(~violated) = 0; clear Ax
            output.constrviolation.lin_ineq{1} = violated;
            output.constrviolation.lin_ineq{2} = violation;
            is_violated = is_violated || any(violated(:));
            max_violation = max(max_violation, max(violation));
            clear violation violated
        end
        if have_lineqconFcn
            Aeqx = Aeq*x_new;
            violated  = abs(Aeqx - beq) > tolCon;
            violation = Aeqx - beq;
            violation(~violated) = 0; clear Aeqx
            output.constrviolation.lin_eq{1} = violated;
            output.constrviolation.lin_eq{2} = violation;
            is_violated = is_violated || any(abs(violated(:)));
            max_violation = max(max_violation, max(violation));
            clear violation violated
        end
        if have_nonlconFcn
            [c, ceq] = feval(nonlcon, x_new);
            output.ConstrfuncCount = output.ConstrfuncCount + 1;
            if ~isempty(ceq)
                violated = abs(ceq) > tolCon;
                ceq(~violated) = 0;
                output.constrviolation.nonlin_eq{1} = violated;
                output.constrviolation.nonlin_eq{2} = ceq;
                is_violated = is_violated || any(violated(:));
                max_violation = max(max_violation, max(abs(ceq)));
                clear violation violated ceq
            end
            if ~isempty(c)
                violated = c > tolCon;
                c(~violated) = 0;
                output.constrviolation.nonlin_ineq{1} = violated;
                output.constrviolation.nonlin_ineq{2} = c;
                is_violated = is_violated || any(violated(:));
                max_violation = max(max_violation, max(c));
                clear violation violated c
            end
            clear c ceq
        end
                
        % adjust output message
        if create_output && exitflag == -3 
            output.message = sprintf(...
                ' No finite function values encountered.\n');
        end
        if ~isfield(output, 'message'), output.message = ''; end % (safeguard)
        if is_violated
            exitflag = -2;
            message = sprintf(...
                [' Unfortunately, the solution is infeasible for the given value ',...
                'of OPTIONS.TolCon of %1.6e\n\t Maximum constraint violation: ',...
                '%1.6e'], tolCon, max_violation);
            clear max_violation
        else
            if exitflag >= 1, message = sprintf('\b\n and'); 
            else message = sprintf('\b\n but');
            end
            message = [message, sprintf([' all constraints are satisfied using ',...
                    'OPTIONS.TolCon of %1.6e.'], tolCon)];
        end
        
        % display or update output structure
        if create_output            
            output.message = sprintf('%s\n%s\n', output.message, message); end
        if do_display && ~do_global_opt
            fprintf(1, '%s\n', message); end
        
        % correct for output possibly wrongfully created above
        if ~create_output, output = []; end
                
    end % finalize  
    
    % optimize globally
    function [sol, fval, exitflag, output] = optimize_globally(varargin)
                
        % first perform error checks
        if isempty(ub) || isempty(lb) || any(isinf(lb)) || any(isinf(ub))
            error('optimize:lbub_undefined',...
                  ['When optimizing globally ([x0] is empty), both [lb] and [ub] ',...
                  'must be non-empty and finite.'])
        end  
        
        % global miniimum (for output function)
        glob_min = inf;
        
        % we can give the popsize in the options structure...
        if isfield(options, 'popsize')
            popsize = options.popsize; 
        % or we use 25*(number of dimensions) individuals by default
        else popsize = 25*numel(lb);
        end
        
        % initialize population of random starting points
        population = repmat(lb,[1,1,popsize]) + ...
            rand(size(lb,1),size(lb,2), popsize).*repmat(ub-lb,[1,1,popsize]);
       
        % get options, and reset maximum allowable function evaluations
        maxiters   = optimget(options, 'MaxIter', 200*numel(lb));
        maxfuneval = optimget(options, 'MaxFunEvals', 1e4);
        MaxFunEval = floor( optimget(options, 'MaxFunEvals', 1e4) / popsize / 1.2);
        options    = optimset(options, 'MaxFunEvals', MaxFunEval);        
        
        % create globalized wrapper for outputfunctions
        have_glob_OutputFcn = false;
        if isfield(options, 'OutputFcn') && ~isempty(options.OutputFcn)
            have_glob_OutputFcn = true;
            glob_OutputFcn = options.OutputFcn;
            options.OutputFcn = @glob_OutputFcn_wrapper;            
        end
        
        % first evaluate output function
        if have_glob_OutputFcn
            optimValues.iteration = 0;
            optimValues.x    = x0;
            optimValues.fval = glob_min;
            optimValues.procedure = 'init';
            optimValues.funcCount = 0;
            if have_nonlconFcn % constrained problems
                optimValues.ConstrfuncCount = 1; end % one evaluation in check_input()
            glob_OutputFcn(x0, optimValues, 'init', varargin{:});
        end
        
        % display header
        if do_display
            if do_extended_display
                fprintf(1, ['  Iter  evals    min f(x)    global min f(x)   max violation\n',....
                    '=====================================================================\n']);
            else
                fprintf(1, ['  Iter  evals    min f(x)    global min f(x)\n',....
                    '============================================\n']);
            end
        end
        
        % slightly loosen options for global method, and
        % kill all display settings
        global_options = optimset(options, ...
            'TolX'   , 1e2 * options.TolX,...
            'TolFun' , 1e2 * options.TolFun,...
            'display', 'off');  
        
        % preserve ConstraintsInObjectiveFunction and TolCon (optimset 
        % deletes these entries)
        if confcn_in_objfcn
            global_options.ConstraintsInObjectiveFunction = confcn_in_objfcn; end        
        if tolCon ~= 1e-6
            options.TolCon = tolCon; end
                
        % initialize loop
        best_fval = inf; iterations = 0; obj_evals = 0; new_x = population(:,:,1);
        sol = NaN(size(lb)); fval = inf; exitflag = []; output = struct; con_evals = 0;        

        % loop through each individual, and use it as initial value
        for ii = 1:popsize
                      
            % optimize current problem
            [sol_i, fval_i, exitflag_i, output_i] = ...
                optimize(funfcn, population(:,:,ii), A, b, Aeq, beq, lb, ub, nonlcon,...
                strictness, global_options, algorithm, varargin{:});
            
            % add number of evaluations and iterations to total
            if ~have_nonlconFcn % unconstrained problems
                obj_evals = obj_evals + output_i.funcCount;
            else % constrained problems
                obj_evals = obj_evals + output_i.ObjfuncCount;
                con_evals = con_evals + output_i.ConstrfuncCount;
            end   
            iterations = iterations + output_i.iterations;
            
            % keep track of the best solution found so far
            if fval_i < best_fval              
                % output values                
                fval = fval_i;   exitflag = exitflag_i;
                sol  = sol_i;    output   = output_i;
                % and store the new best
                best_fval = fval_i;
            end
            
            % reset output structure
            if create_output
                if ~have_nonlconFcn % unconstrained problems
                    output.funcCount = obj_evals;
                else % constrained problems
                    output.ObjfuncCount = obj_evals;
                    output.ConstrfuncCount = con_evals;
                end
                output.iterations = iterations;
            end
            
            % display output so far
            if do_display
                % iter-detailed: include max. constraint violation
                if do_extended_display
                    % do a dummy finalization to get the maximum violation
                    output_j = finalize(sol, output, exitflag_i);
                    max_violation = 0;
                    if have_nonlconFcn
                        max_violation = max([max_violation
                            output_j.constrviolation.nonlin_ineq{2}
                            abs(output_j.constrviolation.nonlin_eq{2})]);
                    end
                    if have_lineqconFcn
                        max_violation = max([max_violation
                            abs(output_j.constrviolation.lin_eq{2})]);
                    end
                    if have_linineqconFcn
                        max_violation = max([max_violation
                            output_j.constrviolation.lin_ineq{2}]);
                    end
                    % print everything
                    fprintf(1, '%4.0d%8.0d%14.4e%15.4e%16.4e\n', ...
                        ii,obj_evals, fval_i,best_fval,max_violation);
                % iter: don't
                else
                    % just print everything
                    fprintf(1, '%4.0d%8.0d%14.4e%15.4e\n', ...
                        ii,obj_evals, fval_i,best_fval);
                end
            end

            % MaxIters & MaxEvals check. The output function may also have 
            % stopped the global optimization
            if (exitflag_i == -1) || (iterations >= maxiters) || ...
               (obj_evals >= maxfuneval)
                % finalize solution                
                [dummy, exitflag] = finalize(sol, output, exitflag);%#ok
                % and break (NOT return; otherwise the last evaluation of 
                % the outputfunction will be skipped)
                break
            end
            
        end % for
   
        % final evaluate output function
        if have_glob_OutputFcn
            optimValues.iteration = iterations;            
            optimValues.procedure = 'optimization complete';  
            optimValues.fval      = fval;
            optimValues.x         = sol;
            optimValues.funcCount = obj_evals;
            if have_nonlconFcn % constrained problems
                optimValues.ConstrfuncCount  = con_evals; end
            glob_OutputFcn(sol, optimValues, 'done', varargin{:});
        end
        
        % check for INF or NaN values. If there are any, finalize 
        % solution and return
        if ~isfinite(fval)    
            [output, exitflag] = finalize(sol, output, -3); return; end
         
        % reset max. number of function evaluations
        options.MaxFunEvals = maxfuneval - obj_evals;
        if have_nonlconFcn % correction for constrained problems
            options.MaxFunEvals = maxfuneval - obj_evals-con_evals; end
        % make sure the display if OFF
        options.Display = 'off';
        % perform the final iteration on the best solution found
        % NOTE: optimize with the stricter options
        [sol, fval, exitflag, output_i] = ...
            optimize(funfcn, sol, A, b, Aeq, beq, lb, ub, nonlcon,...
            strictness, options, algorithm, varargin{:});
        % adjust output
        if create_output
            if ~have_nonlconFcn % unconstrained problems
                output.funcCount = output.funcCount + output_i.funcCount;
            else % constrained problems
                output.ObjfuncCount    = output.ObjfuncCount + output_i.ObjfuncCount;
                output.ConstrfuncCount = output.ConstrfuncCount + output_i.ConstrfuncCount;
            end
            output.iterations = output.iterations + output_i.iterations;
        end
        
        % get the final display right
        if do_display
            fprintf(1, output_i.message); end

        % create temporary message to get the display right
        output.message = output_i.message;
        
        % globalized wrapper for output functions
        function stop = glob_OutputFcn_wrapper(x, optimvalues, state) 
            % only evaluate if the current function value is better than
            % the best thus far found. Also evaluate on on first and last call
            stop = false;
            if (optimvalues.fval <= glob_min) &&...
                    ~any(strcmpi(state, {'done'; 'init'}))
                glob_min = optimvalues.fval;
                stop = glob_OutputFcn(x, optimvalues, state);            
            end
        end
        
    end % optimize_globally
        
end % function

%% Internal Nelder-Mead algorithm

% This algorithm is exactly the same as MATLAB's FMINSEARCH(), except
% that it has several improvements with regard to efficiency; see
%
% Singer & Singer, Applied Numerical Analysis & Computational Mathematics
% Volume 1 Issue 2, Pages 524 - 534,
% "Efficient Implementation of the Nelder-Mead Search Algorithm"
%
% for details. 
function [sol, fval, exitflag, output] = NelderMead(funfcn, x0, options, varargin)

    % Nelder-Mead algorithm control factors
    alpha = 1;     beta  = 0.5;
    gamma = 2;     delta = 0.5;
    a     = 1/20; % (a) is the size of the initial simplex.
    % 1/20 is 5% of the initial values.
    
    % constants
    N                         = numel(x0);
    VSfactor_reflect          = alpha^(1/N);
    VSfactor_expand           = (alpha*gamma)^(1/N);
    VSfactor_inside_contract  = beta^(1/N);
    VSfactor_outside_contract = (alpha*beta)^(1/N);
    VSfactor_shrink           = delta;
    narg                      = nargin;

    % check for output functions
    have_Outputfcn = false;
    if ~isempty(options.OutputFcn)
        have_Outputfcn = true;
        OutputFcn = options.OutputFcn;
    end

    % initial values
    iterations   = -1; sort_inds = (1:N).';
    operation    = 0;  op = 'initial simplex';
    volume_ratio = 1;
    stop = false;

    % parse options
    if (narg == 2) || isempty(options), options = optimset; end
    reltol_x = optimget(options, 'TolX', 1e-4);
    reltol_f = optimget(options, 'TolFun', 1e-4);
    max_evaluations = optimget(options, 'MaxFunEvals', 200*N);
    max_iterations  = optimget(options, 'MaxIter', 1e4);
    display = optimget(options, 'display', 'off');

    % generate initial simplex
    p = a/N/sqrt(2) * (sqrt(N+1) + N - 1) *  eye(N);
    q = a/N/sqrt(2) * (sqrt(N+1) - 1)     * ~eye(N);
    x = x0(:, ones(1,N));
    x = [x0, x + p + q];
   
    % function is known to be ``vectorized''
    f = funfcn(x);    evaluations = N+1;
    
    % first evaluate output function
    if have_Outputfcn
        optimValues.iteration = iterations;
        optimValues.funcCount = evaluations;
        optimValues.fval      = f(1);
        optimValues.procedure = op;
        OutputFcn(x0, optimValues, 'init', varargin{:});
    end

    % sort and re-label initial simplex
    [f, inds] = sort(f);  x = x(:, inds);
    
    % compute initial centroid
    C = sum(x(:, 1:end-1), 2)/N;

    % display header if per-iteration display is selected
    if strcmpi(display, 'iter')
        fprintf(1, ['\n\t\tf(1)\t\tfunc. evals.\toperation\n', ...
            '\t==================================================\n']);
    end

    % main loop
    while true

        % evaluate output function
        if have_Outputfcn
            optimValues.iteration  = iterations;
            optimValues.funcCount  = evaluations;
            optimValues.procedure  = op;
            [optimValues.fval,ind] = min(f);
            stop = OutputFcn(x(:,ind), optimValues, 'iter', varargin{:});
            if stop, break, end
        end

        % increase number of iterations
        iterations = iterations + 1;

        % display string for per-iteration display
        if strcmpi(display, 'iter')
            fprintf(1, '\t%+1.6e', f(1));
            fprintf(1, '\t\t%4.0d', evaluations);
            fprintf(1, '\t\t%s\n', op);
        end
        
        % re-sort function values
        x_replaced = x(:, end);
        if operation == 2  % shrink steps
            [f, inds] = sort(f);   x = x(:, inds);
        else   % non-shrink steps
            inds = f(end) <= f(1:end-1);
            f = [f(sort_inds(~inds)), f(end), f(sort_inds(inds))];
            x = [x(:, sort_inds(~inds)), x_replaced, x(:, sort_inds(inds))];
        end
        
        % update centroid (Singer & Singer are wrong here...
        % shrink & non-shrink steps should be treated the same)
        C = C + (x_replaced -  x(:, end))/N;

        % Algorithm termination conditions
        term_f = abs(f(end) - f(1)) / (abs(f(end)) + abs(f(1)))  < reltol_f;
        fail   = (iterations >= max_iterations) || ...
                 (evaluations >= max_evaluations);
        fail2  = all(~isfinite(f(:))) || all(~isfinite(x(:)));
        term_x = volume_ratio < reltol_x;
        if (term_x || term_f || fail || fail2), break, end

        % non-shrink steps are taken most of the time. Set this as the
        % default operation.
        operation = 1;

        % try to reflect the simplex
        xr = C + alpha*(C - x(:, end));
        fr = funfcn(xr);
        evaluations = evaluations + 1;

        % accept the reflection point
        if fr < f(end-1)
            x(:, end) = xr;
            f(end) = fr;

            % try to expand the simplex
            if fr < f(1)
                xe = C + gamma*(xr - C);
                fe = funfcn(xe);
                evaluations = evaluations + 1;

                % accept expand
                if (fe < f(1))
                    op = 'expand';
                    volume_ratio = VSfactor_expand * volume_ratio;
                    x(:, end) = xe;
                    f(end) = fe;
                    continue;
                end
            end

            % otherwise, just continue
            op = 'reflect';
            volume_ratio = VSfactor_reflect * volume_ratio;
            continue;

        % otherwise, try to contract the simplex
        else

            % outside contraction
            if fr < f(end)
                xc = C + beta*(xr - C);
                insouts = 1;

                % inside contraction
            else
                xc = C + beta*(x(:, end) - C);
                insouts = 2;
            end
            fc = funfcn(xc);
            evaluations = evaluations + 1;

            % accept contraction
            if fc < min(fr, f(end))
                switch insouts
                    case 1
                        op = 'outside contraction';
                        volume_ratio = VSfactor_outside_contract * volume_ratio;
                    case 2
                        op = 'inside contraction';
                        volume_ratio = VSfactor_inside_contract * volume_ratio;
                end
                x(:, end) = xc;
                f(end) = fc;
                continue;

            % everything else has failed - shrink the simplex towards x1
            else                
                % first shrink
                operation = 2;
                xones = x(:, ones(1, N+1));
                x = xones + delta*(x - xones);
                f = funfcn(x);
                evaluations = evaluations + N + 1;
                volume_ratio = VSfactor_shrink * volume_ratio;
                % then evaluate output function
                op = 'shrink';                
            end
            
        end % select next procedure        
    end % main loop

    % evaluate output function
    if have_Outputfcn
        optimValues.iteration = iterations;
        optimValues.funcCount = evaluations;
        optimValues.fval      = f(1);
        optimValues.procedure = 'Final simplex.';
        OutputFcn(x(:,1), optimValues, 'done', varargin{:});
    end

    % final values
    sol  = x(:,1);
    fval = f(1);

    % exitflag
    if (term_x || term_f), exitflag = 1; end % normal convergence
    if stop, exitflag = -1; end              % stopped by outputfunction
    if fail, exitflag = 0; end               % max. iterations or max. func. eval. exceeded
    if fail2, exitflag = -3; end             % everything is non-finite

    % create output structure
    output.iterations = iterations;
    output.funcCount  = evaluations;
    switch exitflag
        case 1
            output.message = sprintf(['Optimization terminated:\n',...
                ' the current x satisfies the termination criteria using OPTIONS.TolX of %d \n',...
                ' and F(X) satisfies the convergence criteria using OPTIONS.TolFun of %d \n'],...
                reltol_x, reltol_f);
        case 0
            if (iterations >= max_iterations)
                output.message = sprintf(['Optimization terminated:\n',...
                    ' Maximum amount of iterations exceeded; \nIncrease ',...
                    '''MaxIters'' option.\n']);
            elseif (evaluations >= max_evaluations)
                output.message = sprintf(['Optimization terminated:\n',...
                    ' Maximum amount of function evaluations exceeded; \nIncrease ',...
                    '''MaxFunevals'' option.\n']);
            end
        case -1
            output.message = ...
                sprintf('Optimization terminated by user-provided output function.\n');
        case -3
            output.message = ...
                sprintf('All function values are non-finite. Exiting...\n');
    end

    % display convergence 
    if ~isempty(options.Display) && ~strcmpi(options.Display, 'off')
        fprintf(1, '\n%s\n', output.message); end

    % make sure the algorithm is correct
    output.algorithm  = 'Nelder-Mead simplex direct search';

end % NelderMead


%% Internal Quasi-Newton L-BFGS algorithm
%FMINLBFGS finds a local minimum of a function of several variables.
%   This optimizer is developed for image registration methods with large
%	amounts of unknown variables.
%
%   Optimization methods supported:
%	- Quasi Newton Broyden�Fletcher�Goldfarb�Shanno (BFGS)
%   - Limited memory BFGS (L-BFGS)
%   - Steepest Gradient Descent optimization.
%
%   [X,FVAL,EXITFLAG,OUTPUT,GRAD] = FMINLBFGS(FUN,X0,OPTIONS)
%
%   Inputs,
%		FUN: Function handle or string which is minimized, returning an
%				error value and optional the error gradient.
%		X0: Initial values of unknowns can be a scalar, vector or matrix
%	 (optional)
%		OPTIONS: Structure with optimizer options, made by a struct or
%				optimset. (optimset doesnot support all input options)
%
%   Outputs,
%		X : The found location (values) which minimize the function.
%		FVAL : The minimum found
%		EXITFLAG : Gives value, which explain why the minimizer stopt
%		OUTPUT : Structure with all important ouput values and parameters
%		GRAD : The gradient at this location
%
%   Extended description of input/ouput variables
%   OPTIONS,
%		OPTIONS.GoalsExactAchieve : If set to 0, a line search method is
%               used which uses a few function calls to do a good line
%               search. When set to 1 a normal line search method with Wolfe
%				conditions is used (default).
%		OPTIONS.GradConstr, Set this variable to true if gradient calls are
%				cpu-expensive (default). If false more gradient calls are
%				used and less function calls.
%	    OPTIONS.HessUpdate : If set to 'bfgs', Broyden�Fletcher�Goldfarb�Shanno
%				optimization is used (default), when the number of unknowns is
%				larger then 3000 the function will switch to Limited memory BFGS,
%				or if you set it to 'lbfgs'. When set to 'steepdesc', steepest
%				decent optimization is used.
%		OPTIONS.StoreN : Number of itterations used to approximate the Hessian,
%			 	in L-BFGS, 20 is default. A lower value may work better with
%				non smooth functions, because than the Hessian is only valid for
%				a specific position. A higher value is recommend with quadratic equations.
%		OPTIONS.GradObj : Set to 'on' if gradient available otherwise finited difference
%				is used.
%     	OPTIONS.Display : Level of display. 'off' displays no output; 'plot' displays
%				all linesearch results in figures. 'iter' displays output at  each
%               iteration; 'final' displays just the final output; 'notify'
%				displays output only if the function does not converge;
%	    OPTIONS.TolX : Termination tolerance on x, default 1e-6.
%	    OPTIONS.TolFun : Termination tolerance on the function value, default 1e-6.
%		OPTIONS.MaxIter : Maximum number of iterations allowed, default 400.
% 		OPTIONS.MaxFunEvals : Maximum number of function evaluations allowed,
%				default 100 times the amount of unknowns.
%		OPTIONS.DiffMaxChange : Maximum stepsize used for finite difference gradients.
%		OPTIONS.DiffMinChange : Minimum stepsize used for finite difference gradients.
%		OPTIONS.OutputFcn : User-defined function that an optimization function calls
%				at each iteration.
%		OPTIONS.rho : Wolfe condition on gradient (c1 on wikipedia), default 0.01.
%		OPTIONS.sigma : Wolfe condition on gradient (c2 on wikipedia), default 0.9.
%		OPTIONS.tau1 : Bracket expansion if stepsize becomes larger, default 3.
%		OPTIONS.tau2 : Left bracket reduction used in section phase,
%		default 0.1.
%		OPTIONS.tau3 : Right bracket reduction used in section phase, default 0.5.
%   FUN,
%		The speed of this optimizer can be improved by also providing
%   	the gradient at X. Write the FUN function as follows
%   	function [f,g]=FUN(X)
%       	f , value calculation at X;
%   	if ( nargout > 1 )
%       	g , gradient calculation at X;
%   	end
%	EXITFLAG,
%		Possible values of EXITFLAG, and the corresponding exit conditions
%		are
%  		1, 'Change in the objective function value was less than the specified tolerance TolFun.';
%  		2, 'Change in x was smaller than the specified tolerance TolX.';
%  		3, 'Magnitude of gradient smaller than the specified tolerance';
%  		4, 'Boundary fminimum reached.';
%  		0, 'Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.';
%  		-1, 'Algorithm was terminated by the output function.';
%  		-2, 'Line search cannot find an acceptable point along the current search';
%
%   Examples
%       options = optimset('GradObj','on');
%       X = fminlbfgs(@myfun,2,options)
%
%   	% where myfun is a MATLAB function such as:
%       function [f,g] = myfun(x)
%       f = sin(x) + 3;
%	    if ( nargout > 1 ), g = cos(x); end
%
%   See also OPTIMSET, FMINSEARCH, FMINBND, FMINCON, FMINUNC, @, INLINE.
%
%   Function written by D.Kroon University of Twente (March 2009)

function [x,fval,exitflag,output,grad]=fminlbfgs(funfcn,x_init,optim)
    
    % Read Optimalisation Parameters
    defaultopt = struct('Display','final','HessUpdate','bfgs','GoalsExactAchieve',1,'GradConstr',true,  ...
        'TolX',1e-6,'TolFun',1e-6,'GradObj','off','MaxIter',400,'MaxFunEvals',100*numel(x_init)-1,  ...
        'DiffMaxChange',1e-6,'DiffMinChange',1e-8,'OutputFcn',[], ...
        'rho',0.0100,'sigma',0.900,'tau1',3,'tau2', 0.1, 'tau3', 0.5,'StoreN',20);
    
    if (~exist('optim','var'))
        optim=defaultopt;
    else
        f = fieldnames(defaultopt);
        for i=1:length(f),
            if (~isfield(optim,f{i})||(isempty(optim.(f{i})))), optim.(f{i})=defaultopt.(f{i}); end
        end
    end
    
    % Initialize the data structure
    data.fval=0;
    data.gradient=0;
    data.fOld=[];
    data.xsizes=size(x_init);
    data.numberOfVariables = numel(x_init);
    data.xInitial = x_init(:);
    data.alpha=1;
    data.xOld=data.xInitial;
    data.iteration=0;
    data.funcCount=0;
    data.gradCount=0;
    data.exitflag=[];
    data.nStored=0;
    data.timeTotal=tic;
    data.timeExtern=0;
    % Switch to L-BFGS in case of more than 3000 unknown variables
    if(optim.HessUpdate(1)=='b')
        if(data.numberOfVariables<3000),
            optim.HessUpdate='bfgs';
        else
            optim.HessUpdate='lbfgs';
        end
    end
    
    if(optim.HessUpdate(1)=='l')
        data.deltaX=zeros(data.numberOfVariables,optim.StoreN);
        data.deltaG=zeros(data.numberOfVariables,optim.StoreN);
        data.saveD=zeros(data.numberOfVariables,optim.StoreN);
    end
    
    exitflag=[];
    
    % Display column headers
    if(strcmp(optim.Display,'iter'))
        disp('     Iteration  Func-count   Grad-count         f(x)         Step-size');
    end
    
    % Calculate the initial error and gradient
    data.initialStepLength=1;
    [data,fval,grad]=gradient_function(data.xInitial,funfcn, data, optim);
    data.gradient=grad;
    data.dir = -data.gradient;
    data.gOld=grad;
    data.fInitial = fval;
    data.fPrimeInitial= data.gradient'*data.dir(:);
    
    
    gNorm = norm(data.gradient,Inf);  % Norm of gradient
    data.initialStepLength = min(1/gNorm,5);
    
    % Show the current iteration
    if(strcmp(optim.Display,'iter'))
        s=sprintf('     %5.0f       %5.0f       %5.0f       %13.6g    ',data.iteration,data.funcCount,data.gradCount,data.fInitial); disp(s);
    end
    
    % Hessian intialization
    if(optim.HessUpdate(1)=='b')
        data.Hessian=eye(data.numberOfVariables);
    end
    
    % Call output function
    % if(call_output_function(data,optim,'init')), exitflag=-1; end
    
    % Start Minimizing
    while(true)
        % Update number of itterations
        data.iteration=data.iteration+1;
        
        % Set current lineSearch parameters
        data.TolFunLnS = eps(max(1,abs(data.fInitial )));
        data.fminimum = data.fInitial - 1e16*(1+abs(data.fInitial));
        
        % Make arrays to store linesearch results
        data.storefx=[]; data.storepx=[]; data.storex=[]; data.storegx=[];
        
        % If option display plot, than start new figure
        if(optim.Display(1)=='p'), figure, hold on; end
        
        % Find a good step size in the direction of the gradient: Linesearch
        if(optim.GoalsExactAchieve==1)
            data=linesearch(funfcn, data,optim);
        else
            data=linesearch_simple(funfcn, data, optim);
        end
        
        % Make linesearch plot
        if(optim.Display(1)=='p');
            plot(data.storex,data.storefx,'r*');
            plot(data.storex,data.storefx,'b');
            
            alpha_test= linspace(min(data.storex(:))/3, max(data.storex(:))*1.3, 10);
            falpha_test=zeros(1,length(alpha_test));
            for i=1:length(alpha_test)
                [data,falpha_test(i)]=gradient_function(data.xInitial(:)+alpha_test(i)*data.dir(:),funfcn, data, optim);
            end
            plot(alpha_test,falpha_test,'g');
            plot(data.alpha,data.f_alpha,'go','MarkerSize',8);
        end
        
        % Check if exitflag is set
        if(~isempty(data.exitflag))
            exitflag=data.exitflag;
            data.xInitial=data.xOld;
            data.fInitial=data.fOld;
            data.gradient=data.gOld;
            break,
        end;
        
        % Update x with the alpha step
        data.xInitial = data.xInitial + data.alpha*data.dir;
        
        % Set the current error and gradient
        data.fInitial =  data.f_alpha;
        data.gradient = data.grad;
        
        % Set initial steplength to 1
        data.initialStepLength = 1;
        
        
        gNorm = norm(data.gradient,Inf);  % Norm of gradient
        
        % Set exit flags
        if(gNorm <optim.TolFun), exitflag=1; end
        if(max(abs(data.xOld-data.xInitial)) <optim.TolX), exitflag=2; end
        if(data.iteration>=optim.MaxIter), exitflag=0; end
        
        % Check if exitflag is set
        if(~isempty(exitflag)), break, end;
        
        % Update the inverse Hessian matrix
        if(optim.HessUpdate(1)~='s')
            % Do the Quasi-Neton Hessian update.
            data = updateQuasiNewtonMatrix_LBFGS(data,optim);
        else
            data.dir = -data.gradient;
        end
        
        % Derivative of direction
        data.fPrimeInitial= data.gradient'*data.dir(:);
        
        % Call output function
        %     if(call_output_function(data,optim,'iter')), exitflag=-1; end
        
        % Show the current iteration
        if(strcmp(optim.Display(1),'i')||strcmp(optim.Display(1),'p'))
            s=sprintf('     %5.0f       %5.0f       %5.0f       %13.6g   %13.6g',data.iteration,data.funcCount,data.gradCount,data.fInitial,data.alpha); disp(s);
        end
        
        % Keep the variables for next iteration
        data.fOld=data.fInitial;
        data.xOld=data.xInitial;
        data.gOld=data.gradient;
    end
    % Set output parameters
    fval=data.fInitial;
    grad=data.gradient;
    x = data.xInitial;
    
    % Reshape x to original shape
    x=reshape(x,data.xsizes);
    
    % Call output function
    % if(call_output_function(data,optim,'done')), exitflag=-1; end
    
    % Make exist output structure
    if(optim.HessUpdate(1)=='b'), output.algorithm='Broyden�Fletcher�Goldfarb�Shanno (BFGS)';
    elseif(optim.HessUpdate(1)=='l'), output.algorithm='limited memory BFGS (L-BFGS)';
    else output.algorithm='Steepest Gradient Descent';
    end
    output.message=getexitmessage(exitflag);
    output.iterations = data.iteration;
    output.funcCount = data.funcCount;
    output.fval = data.fInitial;
    output.stepsize = data.alpha;
    output.directionalderivative = data.fPrimeInitial;
    output.gradient = reshape(data.gradient, data.xsizes);
    output.searchdirection = data.dir;
    output.timeTotal=toc(data.timeTotal);
    output.timeExtern=data.timeExtern;
    oupput.timeIntern=output.timeTotal-output.timeExtern;
    % Display final results
    if(~strcmp(optim.Display,'off'))
        disp('    Optimizer Results')
        disp(['        Algorithm Used: ' output.algorithm]);
        disp(['        Exit message : ' output.message]);
        disp(['        iterations : '  int2str(data.iteration)]);
        disp(['        Function Count : ' int2str(data.funcCount)]);
        disp(['        Minimum found : ' num2str(fval)]);
        disp(['        Intern Time : ' num2str(oupput.timeIntern) ' seconds']);
        disp(['        Total Time : ' num2str(output.timeTotal) ' seconds']);
    end
end

function message=getexitmessage(exitflag)
    switch(exitflag)
        case +1, message='Change in the objective function value was less than the specified tolerance TolFun.';
        case +2, message='Change in x was smaller than the specified tolerance TolX.';
        case +3, message='Magnitude of gradient smaller than the specified tolerance';
        case +4, message='Boundary fminimum reached.';
        case +0, message='Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.';
        case -1, message='Algorithm was terminated by the output function.';
        case -2, message='Line search cannot find an acceptable point along the current search';
        otherwise, message='Undefined exit code';
    end
end

function stop=call_output_function(data,optim,where)
    stop=false;
    if(~isempty(optim.OutputFcn))
        output.iteration = data.iteration;
        output.funcCount = data.funcCount;
        output.fval = data.fInitial;
        output.stepsize = data.alpha;
        output.directionalderivative = data.fPrimeInitial;
        output.gradient = reshape(data.gradient, data.xsizes);
        output.searchdirection = data.dir;
        stop=feval(optim.OutputFcn,...
            reshape(data.xInitial,data.xsizes),output,where);
    end
end

function data=linesearch_simple(funfcn, data, optim)
    % Find a bracket of acceptable points
    data = bracketingPhase_simple(funfcn, data, optim);
    
    if (data.bracket_exitflag  == 2)
        % BracketingPhase found a bracket containing acceptable points;
        % now find acceptable point within bracket
        data = sectioningPhase_simple(funfcn, data, optim);
        data.exitflag = data.section_exitflag;
    else
        % Already acceptable point found or MaxFunEvals reached
        data.exitflag = data.bracket_exitflag;
    end
end

function data = bracketingPhase_simple(funfcn, data,optim)
    % Number of itterations
    itw=0;
    
    % Point with smaller value, initial
    data.beta=0;
    data.f_beta=data.fInitial;
    data.fPrime_beta=data.fPrimeInitial;
    
    % Initial step is equal to alpha of previous step.
    alpha = data.initialStepLength;
    
    % Going up hill
    hill=false;
    
    % Search for brackets
    while(true)
        % Calculate the error registration gradient
        if(optim.GradConstr)
            [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
            fPrime_alpha=nan;
            grad=nan;
        else
            [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
            fPrime_alpha = grad'*data.dir(:);
        end
        
        % Store values linesearch
        data.storefx=[data.storefx f_alpha];
        data.storepx=[data.storepx fPrime_alpha];
        data.storex=[data.storex alpha];
        data.storegx=[data.storegx grad(:)];
        
        % Update step value
        if(data.f_beta<f_alpha),
            % Go to smaller stepsize
            alpha=alpha*optim.tau3;
            
            % Set hill variable
            hill=true;
        else
            % Save current minium point
            data.beta=alpha; data.f_beta=f_alpha; data.fPrime_beta=fPrime_alpha; data.grad=grad;
            if(~hill)
                alpha=alpha*optim.tau1;
            end
        end
        
        % Update number of loop iterations
        itw=itw+1;
        
        if(itw>(log(optim.TolFun)/log(optim.tau3))),
            % No new optium found, linesearch failed.
            data.bracket_exitflag=-2; break;
        end
        
        if(data.beta>0&&hill)
            % Get the brackets around minimum point
            % Pick bracket A from stored trials
            [t,i]=sort(data.storex,'ascend');
            storefx=data.storefx(i);storepx=data.storepx(i); storex=data.storex(i);
            [t,i]=find(storex>data.beta,1);
            if(isempty(i)), [t,i]=find(storex==data.beta,1); end
            alpha=storex(i); f_alpha=storefx(i); fPrime_alpha=storepx(i);
            
            % Pick bracket B from stored trials
            [t,i]=sort(data.storex,'descend');
            storefx=data.storefx(i);storepx=data.storepx(i); storex=data.storex(i);
            [t,i]=find(storex<data.beta,1);
            if(isempty(i)), [t,i]=find(storex==data.beta,1); end
            beta=storex(i); f_beta=storefx(i); fPrime_beta=storepx(i);
            
            % Calculate derivatives if not already calculated
            if(optim.GradConstr)
                gstep=data.initialStepLength/1e6;
                if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
                if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
                [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
                [data,f_beta2]=gradient_function(data.xInitial(:)+(beta+gstep)*data.dir(:),funfcn, data, optim);
                fPrime_alpha=(f_alpha2-f_alpha)/gstep;
                fPrime_beta=(f_beta2-f_beta)/gstep;
            end
            
            % Set the brackets A and B
            data.a=alpha; data.f_a=f_alpha; data.fPrime_a=fPrime_alpha;
            data.b=beta; data.f_b=f_beta; data.fPrime_b=fPrime_beta;
            
            % Finished bracketing phase
            data.bracket_exitflag  = 2; return
        end
        
        % Reached max function evaluations
        if(data.funcCount>=optim.MaxFunEvals), data.bracket_exitflag=0; return; end
    end
end

function data = sectioningPhase_simple(funfcn, data, optim)
    % Get the brackets
    brcktEndpntA=data.a; brcktEndpntB=data.b;
    
    % Calculate minimum between brackets
    [alpha,f_alpha_estimated] = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,data.a,data.b,data.f_a,data.fPrime_a,data.f_b,data.fPrime_b,optim);
    if(isfield(data,'beta')&&(data.f_beta<f_alpha_estimated)), alpha=data.beta; end
    
    
    [t,i]=find(data.storex==alpha,1);
    if((~isempty(i))&&(~isnan(data.storegx(i))))
        f_alpha=data.storefx(i); grad=data.storegx(:,i);
    else
        % Calculate the error and gradient for the next minimizer itteration
        [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
        if(isfield(data,'beta')&&(data.f_beta<f_alpha)),
            alpha=data.beta;
            if((~isempty(i))&&(~isnan(data.storegx(i))))
                f_alpha=data.storefx(i); grad=data.storegx(:,i);
            else
                [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
            end
        end
    end
    
    % Store values linesearch
    data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha];
    
    fPrime_alpha = grad'*data.dir(:);
    data.alpha=alpha;
    data.fPrime_alpha= fPrime_alpha;
    data.f_alpha= f_alpha;
    data.grad=grad;
    
    % Set the exit flag to succes
    data.section_exitflag=[];
    
end

function data=linesearch(funfcn, data, optim)
    
    % Find a bracket of acceptable points
    data = bracketingPhase(funfcn, data,optim);
    
    if (data.bracket_exitflag  == 2)
        % BracketingPhase found a bracket containing acceptable points;
        % now find acceptable point within bracket
        data = sectioningPhase(funfcn, data, optim);
        data.exitflag = data.section_exitflag;
    else
        % Already acceptable point found or MaxFunEvals reached
        data.exitflag = data.bracket_exitflag;
    end
end

function data = sectioningPhase(funfcn, data, optim)
    %
    % sectioningPhase finds an acceptable point alpha within a given bracket [a,b]
    % containing acceptable points. Notice that funcCount counts the total number of
    % function evaluations including those of the bracketing phase.
    
    while(true)
        
        % Pick alpha in reduced bracket
        brcktEndpntA = data.a + min(optim.tau2,optim.sigma)*(data.b - data.a);
        brcktEndpntB = data.b - optim.tau3*(data.b - data.a);
        
        % Find global minimizer in bracket [brcktEndpntA,brcktEndpntB] of 3rd-degree
        % polynomial that interpolates f() and f'() at "a" and at "b".
        alpha = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,data.a,data.b,data.f_a,data.fPrime_a,data.f_b,data.fPrime_b,optim);
        
        % No acceptable point could be found
        if (abs( (alpha - data.a)*data.fPrime_a ) <= data.TolFunLnS), data.section_exitflag = -2; return; end
        
        % Calculate value (and gradient if no extra time cost) of current alpha
        if(~optim.GradConstr)
            [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
            fPrime_alpha = grad'*data.dir(:);
        else
            gstep=data.initialStepLength/1e6;
            if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
            if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
            [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
            [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
            fPrime_alpha=(f_alpha2-f_alpha)/gstep;
        end
        
        % Store values linesearch
        data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha];
        
        % Store current bracket position of A
        aPrev = data.a;
        f_aPrev = data.f_a;
        fPrime_aPrev = data.fPrime_a;
        
        % Update the current brackets
        if ((f_alpha > data.fInitial + alpha*optim.rho*data.fPrimeInitial) || (f_alpha >= data.f_a))
            % Update bracket B to current alpha
            data.b = alpha; data.f_b = f_alpha; data.fPrime_b = fPrime_alpha;
        else
            % Wolfe conditions, if true then acceptable point found
            if (abs(fPrime_alpha) <= -optim.sigma*data.fPrimeInitial),
                if(optim.GradConstr)
                    % Gradient was not yet calculated because of time costs
                    [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
                    fPrime_alpha = grad'*data.dir(:);
                end
                % Store the found alpha values
                data.alpha=alpha; data.fPrime_alpha= fPrime_alpha; data.f_alpha= f_alpha;
                data.grad=grad;
                data.section_exitflag = []; return,
            end
            
            % Update bracket A
            data.a = alpha; data.f_a = f_alpha;  data.fPrime_a = fPrime_alpha;
            
            if (data.b - data.a)*fPrime_alpha >= 0
                % B becomes old bracket A;
                data.b = aPrev; data.f_b = f_aPrev;  data.fPrime_b = fPrime_aPrev;
            end
        end
        
        % No acceptable point could be found
        if (abs(data.b-data.a) < eps), data.section_exitflag = -2; return, end
        
        % maxFunEvals reached
        if(data.funcCount >optim.MaxFunEvals), data.section_exitflag = -1; return, end
    end
end

function data = bracketingPhase(funfcn, data, optim)
    % bracketingPhase finds a bracket [a,b] that contains acceptable points; a bracket
    % is the same as a closed interval, except that a > b is allowed.
    %
    % The outputs f_a and fPrime_a are the values of the function and the derivative
    % evaluated at the bracket endpoint 'a'. Similar notation applies to the endpoint
    % 'b'.
    
    % Parameters of bracket A
    data.a = [];
    data.f_a = [];
    data.fPrime_a = [];
    
    % Parameters of bracket B
    data.b = [];
    data.f_b = [];
    data.fPrime_b = [];
    
    % First trial alpha is user-supplied
    % f_alpha will contain f(alpha) for all trial points alpha
    % fPrime_alpha will contain f'(alpha) for all trial points alpha
    alpha = data.initialStepLength;
    f_alpha = data.fInitial;
    fPrime_alpha = data.fPrimeInitial;
    
    % Set maximum value of alpha (determined by fminimum)
    alphaMax = (data.fminimum - data.fInitial)/(optim.rho*data.fPrimeInitial);
    alphaPrev = 0;
    
    while(true)
        % Evaluate f(alpha) and f'(alpha)
        fPrev = f_alpha;
        fPrimePrev = fPrime_alpha;
        
        % Calculate value (and gradient if no extra time cost) of current alpha
        if strcmpi(optim.GradConstr, 'on')
            [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
            fPrime_alpha = grad'*data.dir(:);
        else
            gstep=data.initialStepLength/1e6;
            if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
            if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
            [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
            [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
            fPrime_alpha=(f_alpha2-f_alpha)/gstep;
        end
        
        % Store values linesearch
        data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha];
        
        % Terminate if f < fminimum
        if (f_alpha <= data.fminimum), data.bracket_exitflag = 4; return; end
        
        % Bracket located - case 1 (Wolfe conditions)
        if (f_alpha > (data.fInitial + alpha*optim.rho*data.fPrimeInitial)) || (f_alpha >= fPrev)
            % Set the bracket values
            data.a = alphaPrev; data.f_a = fPrev;  data.fPrime_a = fPrimePrev;
            data.b = alpha; data.f_b = f_alpha;  data.fPrime_b = fPrime_alpha;
            % Finished bracketing phase
            data.bracket_exitflag  = 2; return
        end
        
        % Acceptable steplength found
        if (abs(fPrime_alpha) <= -optim.sigma*data.fPrimeInitial),
            if(optim.GradConstr)
                % Gradient was not yet calculated because of time costs
                [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
                fPrime_alpha = grad'*data.dir(:);
            end
            % Store the found alpha values
            data.alpha=alpha;
            data.fPrime_alpha= fPrime_alpha; data.f_alpha= f_alpha; data.grad=grad;
            % Finished bracketing phase, and no need to call sectioning phase
            data.bracket_exitflag = [];  return
        end
        
        % Bracket located - case 2
        if (fPrime_alpha >= 0)
            % Set the bracket values
            data.a = alpha; data.f_a = f_alpha;  data.fPrime_a = fPrime_alpha;
            data.b = alphaPrev; data.f_b = fPrev; data.fPrime_b = fPrimePrev;
            % Finished bracketing phase
            data.bracket_exitflag  = 2; return
        end
        
        % Update alpha
        if (2*alpha - alphaPrev < alphaMax )
            brcktEndpntA = 2*alpha-alphaPrev;
            brcktEndpntB = min(alphaMax,alpha+optim.tau1*(alpha-alphaPrev));
            % Find global minimizer in bracket [brcktEndpntA,brcktEndpntB] of 3rd-degree polynomial
            % that interpolates f() and f'() at alphaPrev and at alpha
            alphaNew = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,alphaPrev,alpha,fPrev, ...
                fPrimePrev,f_alpha,fPrime_alpha,optim);
            alphaPrev = alpha;
            alpha = alphaNew;
        else
            alpha = alphaMax;
        end
        
        % maxFunEvals reached
        if(data.funcCount >optim.MaxFunEvals), data.bracket_exitflag = -1; return, end
    end
end

function [alpha,f_alpha]= pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,alpha1,alpha2,f1,fPrime1,f2,fPrime2,optim)
    % finds a global minimizer alpha within the bracket [brcktEndpntA,brcktEndpntB] of the cubic polynomial
    % that interpolates f() and f'() at alpha1 and alpha2. Here f(alpha1) = f1, f'(alpha1) = fPrime1,
    % f(alpha2) = f2, f'(alpha2) = fPrime2.
    
    % determines the coefficients of the cubic polynomial with c(alpha1) = f1,
    % c'(alpha1) = fPrime1, c(alpha2) = f2, c'(alpha2) = fPrime2.
    coeff = [(fPrime1+fPrime2)*(alpha2-alpha1)-2*(f2-f1) ...
        3*(f2-f1)-(2*fPrime1+fPrime2)*(alpha2-alpha1) (alpha2-alpha1)*fPrime1 f1];
    
    % Convert bounds to the z-space
    lowerBound = (brcktEndpntA - alpha1)/(alpha2 - alpha1);
    upperBound = (brcktEndpntB - alpha1)/(alpha2 - alpha1);
    
    % Swap if lowerbound is higher than the upperbound
    if (lowerBound  > upperBound), t=upperBound; upperBound=lowerBound; lowerBound=t; end
    
    % Find minima and maxima from the roots of the derivative of the polynomial.
    sPoints = roots([3*coeff(1) 2*coeff(2) coeff(3)]);
    
    % Remove imaginaire and points outside range
    sPoints(imag(sPoints)~=0)=[];
    sPoints(sPoints<lowerBound)=[]; sPoints(sPoints>upperBound)=[];
    
    % Make vector with all possible solutions
    sPoints=[lowerBound sPoints(:)' upperBound];
    
    % Select the global minimum point
    [f_alpha,index]=min(polyval(coeff,sPoints)); z=sPoints(index);
    
    % Add the offset and scale back from [0..1] to the alpha domain
    alpha = alpha1 + z*(alpha2 - alpha1);
    
    % Show polynomial search
    if(optim.Display(1)=='p');
        vPoints=polyval(coeff,sPoints);
        plot(sPoints*(alpha2 - alpha1)+alpha1,vPoints,'co');
        plot([sPoints(1) sPoints(end)]*(alpha2 - alpha1)+alpha1,[vPoints(1) vPoints(end)],'c*');
        xPoints=linspace(lowerBound/3, upperBound*1.3, 50);
        vPoints=polyval(coeff,xPoints);
        plot(xPoints*(alpha2 - alpha1)+alpha1,vPoints,'c');
    end
    
end

function [data,fval,grad]=gradient_function(x,funfcn, data, optim)
    % Call the error function for error (and gradient)
    if ( nargout <3 )
        timem=tic;
        fval=funfcn(reshape(x,data.xsizes));
        data.timeExtern=data.timeExtern+toc(timem);
        data.funcCount=data.funcCount+1;
    else
        if(strcmp(optim.GradObj,'on'))
            timem=tic;
            [fval, grad]=feval(funfcn,reshape(x,data.xsizes));
            data.timeExtern=data.timeExtern+toc(timem);
            data.funcCount=data.funcCount+1;
            data.gradCount=data.gradCount+1;
        else
            % Calculate gradient with forward difference if not provided by the function
            grad=zeros(length(x),1);
            fval=funfcn(reshape(x,data.xsizes));
            gstep=data.initialStepLength/1e6;
            if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
            if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
            for i=1:length(x),
                x_temp=x; x_temp(i)=x_temp(i)+gstep;
                timem=tic;
                [fval_g]=feval(funfcn,reshape(x_temp,data.xsizes)); data.funcCount=data.funcCount+1;
                data.timeExtern=data.timeExtern+toc(timem);
                grad(i)=(fval_g-fval)/gstep;
            end
        end
        grad=grad(:);
    end
end

function data = updateQuasiNewtonMatrix_LBFGS(data,optim)
    % updates the quasi-Newton matrix that approximates the inverse to the Hessian.
    % Two methods are support BFGS and L-BFGS, in L-BFGS the hessian is not
    % constructed or stored.
    % Calculate position, and gradient diference between the
    % itterations
    deltaX=data.alpha* data.dir;
    deltaG=data.gradient-data.gOld;
    
    if ((deltaX'*deltaG) >= sqrt(eps)*max( eps,norm(deltaX)*norm(deltaG) ))
        
        if(optim.HessUpdate(1)=='b')
            % Default BFGS as described by Nocedal
            p_k = 1 / (deltaG'*deltaX);
            Vk = eye(data.numberOfVariables) - p_k*deltaG*deltaX';
            % Set Hessian
            data.Hessian = Vk'*data.Hessian *Vk + p_k * (deltaX*deltaX.');
            % Set new Direction
            data.dir = -data.Hessian*data.gradient;
        else
            % L-BFGS with scaling as described by Nocedal
            
            % Update a list with the history of deltaX and deltaG
            data.deltaX(:,2:optim.StoreN)=data.deltaX(:,1:optim.StoreN-1); data.deltaX(:,1)=deltaX;
            data.deltaG(:,2:optim.StoreN)=data.deltaG(:,1:optim.StoreN-1); data.deltaG(:,1)=deltaG;
            
            data.nStored=data.nStored+1; if(data.nStored>optim.StoreN), data.nStored=optim.StoreN; end
            
            % Initialize variables
            a=zeros(1,data.nStored);
            p=zeros(1,data.nStored);
            
            q = data.gradient;
            for i=1:data.nStored
                p(i)= 1 / (data.deltaG(:,i)'*data.deltaX(:,i));
                a(i) = p(i)* data.deltaX(:,i)' * q;
                q = q - a(i) * data.deltaG(:,i);
            end
            % Scaling of initial Hessian (identity matrix)
            p_k = data.deltaG(:,1)'*data.deltaX(:,1) / sum(data.deltaG(:,1).^2);
            
            % Make r = - Hessian * gradient
            r = p_k * q;
            for i=data.nStored:-1:1,
                b = p(i) * data.deltaG(:,i)' * r;
                r = r + data.deltaX(:,i)*(a(i)-b);
            end
            
            % Set new direction
            data.dir = -r;
        end
    end
end
