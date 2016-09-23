function [V1, V2, params, exitflag, output] = central_body_flyby(r1vec, r2vec, tf, m, ...
        rpmin, confcn_r1, confcn_r2, muC, scope)
% CENTRAL_BODY_FLYBY  

% exitflag 
%    1      problem was sucessfully solved
%   -1      the problem has no solution
%   -2      the algorithm failed to find a feasible solution;
%           returned solution will be unreliable
%   -3      the algorithm failed to find a good initial value.
    
    % last edited 14/Nov/2009
         
    %% INITIALIZE
    
    % NOTE: the input argument [m] is currently not used, simply because I
    % don't know how to solve that problem. This is a possible future
    % extention ^_^
    
    % parse input
    r1       = sqrt(r1vec*r1vec.');             % magnitude of r1
    r2       = sqrt(r2vec*r2vec.');             % magnitude of r2   
    dotprod  = r1vec*r2vec.';                   % dot product r1·r2
    dth      = real(acos(dotprod/r1/r2));       % total turn angle
    rpmax    = min(r1,r2);                      % maximum allowable [rp]    
    longway  = sign(tf); tf = abs(tf);          % solve for the longway if [tf] < 0
    if (longway < 0), dth = 2*pi - dth; end     % otherwise, solve the shortway
    r1unit   = r1vec/r1;                        % unit vector of r1vec        
    r2unit   = r2vec/r2;                        % unit vector of r2vec
    crsprod  = cross(r1vec, r2vec, 2);          % normal vector to both r1vec and r2vec
    mcrsprd  = sqrt(crsprod*crsprod.');         % magnitude of that normal
    th1unit  = cross(crsprod/mcrsprd, r1unit);  % unit vectors in the tangential-directions
    th2unit  = cross(crsprod/mcrsprd, r2unit);    
    maxiters = 5;                               % maximum allowable iterations to find [x0]
    
    % initial output is pessimistic
    V1 = NaN(1, 3);   V2 = V1;  exitflag = -1;
    params.DeltaV = NaN;    params.theta1 = NaN;   
    params.dth    = dth;    params.tf1    = NaN;    
    output = struct(...
        'ObjfuncCount'   , 0,...
        'ConstrfuncCount', 0,...
        'iterations'     , 0,...
        'algorithm'      , [],...
        'constrviolation', struct('nonlin_eq', {false, 0}),...
        'message'        , 'Problem has no solution.');
    
    %% UPPER AND LOWER LIMITS (ALSO -- FIRST CHECK ON SOLVABILITY)
        
    % define upper- and lower limits of the search space
    % To that end, first check whether the maximum [rpmax] allows for
    % solutions. If not find a new [rpmax] that *will* allow solutions. 
    lower_lim = acos(rpmin/r1);       % follows from [e1 >= 0]
    upper_lim = dth - acos(rpmin/r2); % follows from [e2 >= 0]
    if lower_lim > upper_lim
        % this equation follows from solving 
        %   acos(rpmin/r1) = dth - acos(rpmin/r2) 
        rpmin = r2*r1*sin(dth) / sqrt(r1^2 + r2^2 - 2*r1*r2*cos(dth));
        % slightly offset if to get at least some solutions
        rpmin = rpmin*1.01; % plus 1%
        % re-compute the limits
        lower_lim = acos(rpmin/r1);
        upper_lim = dth - acos(rpmin/r2);
    end
    % naturally, this new limit may fall above the maximum allowable limit,
    % in which case there are no solutions
    if (rpmin > rpmax), return, end
    % also, there may be no solution to this equation. Also in that case
    % there can be no solutions
    if (lower_lim > upper_lim), return, end
    
    % otherwise, the lower- and upper limits are simply:
    LB = [  0, rpmin];
    UB = [dth, rpmax];
        
    %% FIND INITIAL ESTIMATE (ALSO -- SECOND CHECK ON SOLVABILITY)
    
    % start with the minimum allowable [rp]. Initialize the loop
    rp = rpmin;             iteration = 0;
    tf_too_short = false;   tf_too_long = false;
    
rp = rpmax-0.1*150e6;
    
    % Several of the limits must be slightly offset when evaluating the
    % time-of-flight function; when it is evaluated at the exact end-points, 
    % it may be that the function value or its derivative is [inf].
    pb_offset  = 1e-4; % offset in case of near-parabolic limits
    nrm_offset = 1e-6; % offset for all other cases
    
    % slightly offset the limits (AT these limits, the constraint e < 0 is
    % still slightly violated)
    lower_lim = lower_lim + nrm_offset;
    upper_lim = upper_lim - nrm_offset;    
    
    % loop until a good initial value has been found        
    while true
        
        % increase iterations
        iteration = iteration + 1;
        
        % Something else may go wrong
        if iteration > maxiters
            exitflag = -3;  % failure
            output.message = 'Failed to find an initial value.';
            return
        end
        
        % these points are singularities
        singul1 = acos(2*rp/r1-1);      singul3 = acos(2*rp/r2-1);
        singul2 = 2*pi - singul1;       singul4 = 2*pi - singul3;
        parab1  = [singul1,singul2];    parab2  = [singul3,singul4];
        
        % check if any are larger than pi, and fall in the current interval
        parab1 = parab1(parab1 > pi & parab1 > lower_lim & parab1 < upper_lim);
        parab2 = dth - parab2(parab2 > pi & parab2 > lower_lim & parab2 < upper_lim);        
        % sort them
        if (parab1 > parab2), t = parab1; parab1 = parab2; parab2 = t; end
        
        % (re-)define time of flight function, its derivative w.r.t. [rp],
        % and the Delta-V equation (re-define, because we have a different
        % [rp] at each iteration)
        tof_eq     = @(th) tf - TOF(th,rp,r1,false,muC) - TOF(dth-th,rp,r2,false,muC);
        dtof_eqdrp = @(th) compute_dtfdrp(th, rp);
        dtof_eqdth = @(th) compute_dtfdth(th, rp);
        DV_eq      = @(th) DeltaV_equation([th, rp]);
        
        % split the time of flight equation into 1, 2, or 3 parts,
        % depending on the presence of singularities in the interval      
        lower_lims = [lower_lim; parab1.'+pb_offset; parab2.'+pb_offset];
        upper_lims = [parab1.'-pb_offset; parab2.'-pb_offset; upper_lim];
        limits     = [lower_lims, upper_lims];
                        
        % evaluate the function at the terminal and center points        
        f_lower_end    = arrayfun(tof_eq, lower_lims);
        f_upper_end    = arrayfun(tof_eq, upper_lims);
        f_center_point = arrayfun(tof_eq, (upper_lims + lower_lims)/2);
        fun_vals       = [f_lower_end, f_center_point, f_upper_end];
        
        % check if any of the three function values has opposite sign to
        % the other two
        good_intervals = abs(sign(f_lower_end) + sign(f_center_point) + sign(f_upper_end)) == 1;
        good_limits    = limits(good_intervals, :);
        
        % For the first iteration, check if the TOF-equation lies entirely 
        % below zero for the smallest allowable value for [rp]. If that is 
        % the case, the problem probably has no solution, and we can return
        if iteration == 1            
            % they ALL lie below zero
            if all(fun_vals(:) < 0)                
                % there can still be roots; get the roots of 
                % [dtf/dth] to check for maxima      
                extrema = [];
                for ii = 1:numel(f_lower_end)
                    % get the roots (avoid the exact end-points)
                    deriv_roots = FindRealRoots(dtof_eqdth,...
                        limits(ii,1)+pb_offset, limits(ii,2)-pb_offset, 10);
                    % evaluate the function at these roots                    
                    if ~isempty(deriv_roots), extrema = arrayfun(tof_eq, deriv_roots); end
                    % if one of these maxima lies above zero, we know that
                    % there is a root in this interval. In that case, just 
                    % skip to finding the roots
                    if any(extrema > 0)
                        good_limits = limits(ii, :);
                        good_intervals = true;
                        % exit the loop
                        break;
                        
                    % However, if all of these maxima also lie below zero,
                    % there is definitely no solution -- we can return.
                    else
                        exitflag = -1; % problem infeasible
                        output.message = 'The given travel time is too short.';
                        return;
                    end
                end
            end            
        end
        
        % If there are good intervals, we know there is at least one 
        % root in one of these intervals. Find all roots of this 
        % interval, and use the best one as the initial value
        if any(good_intervals)            
            % get these roots (avoid the exact end-points)
            roots = FindRealRoots(tof_eq,good_limits(1,1)+pb_offset,good_limits(1,2)-pb_offset,10);
            % finding roots may fail - return in that case
            if isempty(roots), exitflag = -3; return; end
            % evaluate the Delta-V equation at each of these roots
            DeltaV = arrayfun(DV_eq, roots);
            % the one that yields the lowest Delta-V is the initial value
            [minDV, index] = min(DeltaV); %#ok          
            x0 = [roots(index), rp];
            % finding the initial value is complete
            break
           
        % otherwise, we have to shift [rp] and do it all over again
        else            
            % NOTE: the minimum function value must be POSITIVE during
            % the first iteration; otherwise, the next [rp] will always
            % fall below the minimum allowable [rpmin].
            if (iteration == 1), fun_vals(fun_vals < 0) = inf; end
            % Find the minimum function value, and the proper interval
            [minimum_fval_x, index_x] = min(abs(fun_vals),[],1);
            [minimum_fval  , index_y] = min(abs(minimum_fval_x),[],2);
            % the sign is important for consecutive iterations
            minimum_fval = minimum_fval * sign(fun_vals(index_x(index_y),index_y));
            min_limits = limits(index_x(index_y), :);      
            % use the center point point for derivative w.r.t. [rp]
            % (this point is safer to use than the end-points, since
            % the values of the derivatives there can attain
            % extremely large values)
            dtfdrp = dtof_eqdrp(sum(min_limits)/2);
            % do a single "educated-guess" Newton-Raphson step 
            rp = rp - minimum_fval/dtfdrp;
            % [rp] must remain larger than [rpmin] and smaller than [rpmax]
            if rp > rpmax
                % first time round, just give it a nudge
                rp = rpmax - (rpmax-rpmin)/10;
                % the algorithm may try to increase [rp] AGAIN, in which case
                % we know the given travel time is simply too long
                if tf_too_long
                    exitflag = -1; % problem infeasible
                    output.message = 'The given travel time is too long.';
                    return
                end
                % save its state
                tf_too_long = true;    
            end
            if rp < rpmin
                % first time round, just give it a nudge
                rp = rpmin + (rpmax-rpmin)/10;                
                % the algorithm may try to decrease [rp] AGAIN, in which case
                % we know the given travel time is simply too short
                if tf_too_short
                    exitflag = -1; % problem infeasible
                    output.message = 'The given travel time is too short.';
                    return
                end
                % save its state
                tf_too_short = true;                
            end
        end
            
        % adjust the limits for the next iteration
        lower_lim = acos(rp/r1) + nrm_offset;
        upper_lim = dth - acos(rp/r2) - nrm_offset;
                
    end
    
    % small wrapper function to return only [dtfdrp]
    function dtfdrp = compute_dtfdrp(theta, rpp)
        % evaluate the function 
        [dontcare, grad_tf1] = TOF(theta, rpp, r1, true, muC);%#ok
        [dontcare, grad_tf2] = TOF(dth-theta, rpp, r2, true, muC);%#ok
        % return the derivative
        dtfdrp = -grad_tf1(2)-grad_tf2(2);
    end
    
    % small wrapper function to return only [dtfdrp]
    function dtfdth = compute_dtfdth(theta, rpp)
        % evaluate the function 
        [dontcare, grad_tf1] = TOF(theta, rpp, r1, true, muC);%#ok
        [dontcare, grad_tf2] = TOF(dth-theta, rpp, r2, true, muC);%#ok
        % return the derivative
        dtfdth = -grad_tf1(1)+grad_tf2(1);
    end
    
    %% OPTIMIZE THE PROBLEM
        
    % select scope
    if strcmpi(scope, 'global')
        maxiters = 5;    % global problems need to be FAST
    elseif strcmpi(scope, 'local')
        maxiters = 250;  % local problems need to be ACCURATE
    end
    
    % set options accordingly
    % NOTE: don't use optimset (speed)
    options = struct(...
        'Display'    , 'off', ...   % make sure it's *OFF*
        'TolCon'     , 1e-3,...     % about 10 minutes tolerance on the TOF-constraint
        'TolFun'     , 1e-4,...     % compute DeltaV to within 0.1 m/s
        'TolX'       , 1e-7, ...    % needs to be lower due to [sin]-coordinate transformation
        'MaxFunEvals', 250,...      % included for safety
        'MaxIter'    , maxiters,... % ditto
        'Algorithm'  , 'active-set');... % slightly faster, but less robust    
    options.ConstraintsInObjectiveFunction = 5; % evaluate contraints inside the objective function
        
% FMINCON() is much faster and more robust for the CBF-problem than
% OPTIMIZE(). However, FMINCON() gets stuck in an infinite loop inside
% NLCONST() (the loop from line 609ish) sometimes (because of 
% infinities in DV), so until that problem gets fixed it's simply not 
% very safe to use.

    try
        % try optimizing with FMINCON()  
        % NOTE: for global problems, active-set usually finds feasible solutions
        % in less iterations. It is far less robust though...        
        [solution, fval,ig, output] = ...
                fmincon(@DeltaV_equation, x0, [],[],[],[], LB,UB, @nonlcon, options);            
                     
%         % local optimizations are allowed far more iterations. Also, both
%         % interior-point and active-set are used, to see which one yields
%         % the minimum constraint violation
%         if  strcmpi(scope, 'local')
%             options.Algorithm = 'interior-point';
%             % optimize again with FMINCON()
%              [solution2, fval2,ig, output2] = ...
%                 fmincon(@DeltaV_equation, x0, [],[],[],[], LB,UB, @nonlcon, options);
%             % select the better solution of the two
%             if (output2.constrviolation < output.constrviolation) ||...
%                (output2.constrviolation == output.constrviolation && fval2 < fval)
%                 solution = solution2; output = output2;            
%             end            
%         end  
        
    catch %#ok        
        % optimize using constrained FMINSEARCH() via OPTIMIZE() if
        % FMINCON() fails
        [solution, dont,care, output] = ...
            optimize(@DeltaV_equation, x0, [],[],[],[], LB,UB, [],[], options, 'fminsearch');%#ok
    end
    
    % objective/constraint function: Delta-V
    function [DV, V1, V2, tf1, c, ceq] = DeltaV_equation(X,varargin)        
        % split input
        theta = X(1); rp = X(2);
        % limits on [theta] and [rp]
        lim1 = acos(rp/r1);   lim2 = dth - acos(rp/r2);        
        % calculate semi-major axes
        [tf1, dontcare, a1] = TOF(theta, rp, r1, false, muC);%#ok
        [tf2, dontcare, a2] = TOF(dth - theta, rp, r2, false, muC);%#ok
        % calculate Delta-V (at pericenter only)
        DeltaV1 = sqrt(2/rp - 1/a1);  DeltaV2 = sqrt(2/rp - 1/a2);        
        DV      = sqrt(muC) * abs(DeltaV2 - DeltaV1);             
        % compute terminal speeds
        mag_V1 = sqrt(muC*(2/r1 - 1/a1)); % NOTE: *NOT* equal to DeltaV1,2!!
        mag_V2 = sqrt(muC*(2/r2 - 1/a2));
        
        % calculate eccentricities
        e1 = 1 - rp/a1;   e2 = 1 - rp/a2;
        % flight-path angles
        gamma1 = atan2(e1*sin(theta), 1+e1*cos(theta));
        gamma2 = atan2(e2*sin(dth-theta), 1+e2*cos(dth-theta));         
        % compute terminal velocities               
        % (the minus signs still go wrong sometimes..)
        Vr1 = -mag_V1*r1unit*sin(gamma1);   Vtheta1 = longway * mag_V1*th1unit*cos(gamma1);
        Vr2 = +mag_V2*r2unit*sin(gamma2);   Vtheta2 = longway * mag_V2*th2unit*cos(gamma2);
        V1  = Vr1 + Vtheta1;                V2      = Vr2 + Vtheta2;   
        % evaluate boundary functions
        c_1 = confcn_r1(V1);  c_1 = c_1(:);
        c_2 = confcn_r2(V2);  c_2 = c_2(:);
        
        % nonlinear equality constraints
        ceq = (tf - tf1 - tf2)*1e2; % for physical solutions (with some scaling)               
        % nonlinear inequality constraints
        % NOTE: [rp]-constraints are already included via LB, UB
        if strcmpi(scope, 'global')
        c = [ lim1 - theta           % [theta] must always be larger than [lim1]
             theta - lim2            % and smaller than [lim2] 
              sum(c_1(c_1 > 1e-3))   % constraints at first body
              sum(c_2(c_2 > 1e-3))]; % constraints at second body
          % NOTE: taking the sum reduces the required amount of function
          % evaluations in FMINCON to get the gradients of the constraints.
          % Less accurate, but way faster
        elseif strcmpi(scope, 'local')
            c = [ lim1 - theta  % [theta] must always be larger than [lim1]
                theta - lim2    % and smaller than [lim2]
                c_1             % constraints at first body
                c_2];           % constraints at second body
        end       
    end % DeltaV & constraint equation
        
    % wrapper function to het the constraints into FMINCON()
    function [c,ceq] = nonlcon(X,varargin)
         [ig,ig,ig,ig, c, ceq] = DeltaV_equation(X);%#ok
    end
    
    %% OUTPUT SOLUTION
    
    % rename solution
    theta1 = solution(1); % theta1 is a an output argument
    rp     = solution(2);
    
    % get DeltaV and the accompanying semi-major axes and eccentricities
    [DeltaV, V1, V2, tf1] = DeltaV_equation(solution);
    
    % collect output arguments     
    params.DeltaV = DeltaV;        % deltaV used at pericenter
    params.rp     = rp;            % pericenter distance
    rp1_vec = progress_orbit(tf1, [r1vec,V1],muC); % actual statevector to point of closest approach         
    rp2_vec = progress_orbit(tf1-tf, [r2vec,V2],muC); % also from second point
    params.rp_vec = rp1_vec(1:3);  % only position
    params.Vp_vec = [rp1_vec(4:6); % velocity from first body's orbit
                     rp2_vec(4:6)];% velocity from last body's orbit
    params.dth    = dth;           % total turn angle
    params.tf1    = tf1;           % time of flight from first body to pericenter
    params.tf2    = tf - tf1;      % time of flight from pericenter to last body
    params.theta2 = dth-theta1;    % turn angle for first trajectory (always positive)
    params.theta1 = -theta1;       % turn angle for first trajectory (always negative)
    params.constrviolation = 0;    % output the constraint violation
        
    % set appropriate exitflag
    % ([output] is the same one as the one from OPTIMIZE() or FMINCON())
    exitflag = 1;
    if isfield(output, 'constrviolation')
        % OPTIMIZE()
        if isfield(output.constrviolation, 'nonlin_eq') && ...
                output.constrviolation.nonlin_eq{2} > options.TolCon
            exitflag = -2;
            params.constrviolation = output.constrviolation.nonlin_eq{2};
        % FMINCON()
        elseif output.constrviolation > options.TolCon
            exitflag = -2;
            params.constrviolation = output.constrviolation;
        end
    end
    
end % Central body flyby

%% TIME-OF-FLIGHT SUBFUNCTION

% This subfunction has also been compiled in an attempt to increase the
% execution time of the central body flyby targeter. If you wish to use the
% M-version below, simply remove the "_2" from its name.

% this function determines the semi-major axis, eccentricity and time
% of flight for given values of [theta], [r] and [rp].    
function [tf, grad_tf, a, grad_a, e] = TOF_2(theta, rp, r, gradients, muC)

    % set gradients as empty when they are NOT required
    if ~gradients, grad_tf = []; grad_a = []; end
        
    % default output
    tf = inf; 

    % precompute some quantities
      costh = cos(theta);            sinth = sin(theta);
    costhp1 = 1 + costh;              rpsq = rp*rp;
        rsq = r*r;                     rrp = r*rp;
      denom = r*costhp1 - 2*rp;     if denom == 0, denom = realmin; end

    % compute semi-major axis
    a = (rrp*costh - rpsq)/denom; if a == 0, a = realmin; end

    % and its gradient
    % (checked with central differences)
    if gradients
        dadth  = sinth*(r*rpsq-rsq*rp)/denom^2;
        dadrp  = (costhp1*(rsq*costh-2*rrp)+2*rpsq)/denom^2;
        grad_a = [dadth; dadrp];
    end

     % compute eccentricity
    e = 1 - rp/a;

    % pre-compute some more quantities
    sqrt1meo1pe = sqrt( abs(1-e)/(1+e) );
            asq = a*a;
          ep1sq = (e+1)^2;

    % gradient of eccentricity
    % (checked with central differences)
    if gradients
        dedth = rp*dadth/asq;
        dedrp = (rp*dadrp - a)/asq;
    end
    
    % split the different cases        
    elliptic   = e < 1;
    parabolic  = e == 1;
    hyperbolic = e > 1;  
    
    % check if e > 0, rp > 0 and/or theta < acos(-1/e). If these conditions
    % are not true, set everything to [inf] and return
    if (e < 0) || ~isfinite(theta) || ~isfinite(rp) || (rp < 0) || ...
       (hyperbolic && abs(theta) > acos(-1/e))
        tf = inf; grad_tf = [inf; inf]; e = inf;
        a = inf;  grad_a = [inf; inf]; return
    end
    
    % elliptic case
    % (derivatives all checked with central differences)
    if elliptic
        % compute mean motion 
        n       = sqrt(muC/a^3);
        % compute the time of flight
        E       = 2*atan2(sqrt1meo1pe*sin(theta/2), cos(theta/2));
        sinE    = sin(E);   cosE = cos(E);   tantho2 = tan(theta/2);            
        tf      = (E - e*sin(E))/n;
        % and its gradient
        if gradients
            dndth   = -3/2*sqrt(muC/a^5)*dadth;
            dndrp   = dndth*dadrp/dadth;
            dEdth   = 2/(1 + (sqrt1meo1pe*tantho2)^2)*...;
                (sqrt1meo1pe/2/cos(theta/2)^2-tantho2*dedth/ep1sq/sqrt1meo1pe);
            dEdrp   = -2*tantho2*dedrp/(1 + (sqrt1meo1pe*tantho2)^2)/sqrt1meo1pe/ep1sq;
            dtfdth  = (dEdth*(1 - e*cosE) - dedth*sinE - tf*dndth)/n;
            dtfdrp  = (dEdrp*(1 - e*cosE) - dedrp*sinE - tf*dndrp)/n;
            grad_tf = [dtfdth; dtfdrp];
            % scale to days
            grad_tf = grad_tf /86400;
        end
        % scale to days
        tf = tf / 86400;
        
    % parabolic case
    % (derivatives all checked with central differences)
    elseif parabolic
        % compute mean motion        
        n = sqrt(muC/8/rp^3);
        % compute time of flight
        tantho2 = tan(theta/2);
        M = (tantho2 + (tantho2^3)/3)/2;
        tf = M/n;
        % and its gradient
        if gradients
            dndrp  = -3*muC/16/sqrt(rp^5*muC/8);
            dMdth  = 1/4/cos(theta/2)^2 * (1 + tantho2^2);
            dtfdrp = -M*dndrp/n/n;
            dtfdth = dMdth/n;
            % scale to days
            grad_tf = [dtfdth; dtfdrp]/86400;
        end
        % scale to days
        tf = tf / 86400;

     % hyperbolic case
     % (derivatives all checked with central differences)
     elseif hyperbolic
         % compute mean motion 
        n  = sqrt(muC/(-a)^3);
        % compute the time of flight
        tantho2 = tan(theta/2);
        F  = 2*atanh( sqrt1meo1pe*tantho2 );
        tf = (e*sinh(F) - F)/n;            
        % and its gradient
        if gradients
            dndth   = +3/2*sqrt(muC/(-a)^5)*dadth;
            dndrp   = +3/2*sqrt(muC/(-a)^5)*dadrp;
            dFdth   = 1/(1-(sqrt1meo1pe*tantho2)^2) * ...
                (tantho2*2*dedth/ep1sq/sqrt1meo1pe + sqrt1meo1pe*(1 + tantho2^2));
            dFdrp   = tantho2*2*dedrp/ep1sq/sqrt1meo1pe/(1-(sqrt1meo1pe*tantho2)^2);
            dtfdth  = (dedth*sinh(F) + dFdth*(e*cosh(F)-1) - tf*dndth)/n;
            dtfdrp  = (dedrp*sinh(F) + dFdrp*(e*cosh(F)-1) - tf*dndrp)/n;
            grad_tf = [dtfdth; dtfdrp];
            % scale to days
            grad_tf = grad_tf/86400;
        end
        % scale to days
        tf = tf / 86400;
    end

end % TOF
