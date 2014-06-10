function [reachable, encounter_times, min_dist, rel_speed] = ...
        minimum_distance_exposins(exposin, R, secs, t0, tend, muC, params)
    
    %% INITIALIZE

    % initialize output values (it's easiest to be pessimistic)
    reachable       = false(size(secs,1),1);
    encounter_times = NaN(size(secs,1),1);
    min_dist        = inf(size(secs,1),1);
    rel_speed       = min_dist;
    
    % parameters & defaults
    threshold = 0.01 * 1.495978706910000e+008; % default threshold; 0.01 AU    
    all_inds  = 1:size(secs,1);% vector containing all the indices (avoids FIND())    
    
    % errortrap (ALL arguments must be given!)
    error(nargchk(6,7,nargin));%#ok
    
    % times are given in days; convert to seconds
    % (everything else is given in seconds)
    t0   = t0 * 86400;
    tend = tend * 86400;
    
    % total time window
    timespan = tend - t0;
    
    %% MAKE CONSTANT-TIMESTEP CUBIC SPLINES INTERPOLATION OF EXPOSIN
    
    % CONTENTS OF VECTOR [exposin] :
    %
    %   exposin =  [k0,k1,k2,phi,tf,N, dth,gamma1,gamma_m,gamma_M]
    %
    % (See also lambert_low_exposins.m)   
    
    % extract relevant coefficients from ExpoSin
    [k0, k1, k2, phi, dth] = deal(exposin(1),exposin(2),exposin(3),exposin(4),exposin(7));

    % get some positions of the exposin in the time window. We need high
    % accuracy, because eventually we would like constant TIME steps; not
    % steps in [theta]. We have to use steps in [theta] here, as computing 
    % steps in [t] would require a root-finding procedure, with a quadrature 
    % evaluation *at each function evaluation*; pretty costly, and pretty 
    % difficult to implement correctly. The following method is just easier:

    % initialize loop
    f_try = 0:0.01:exposin(7);              t_try = zeros(size(f_try));   
    coordsp_try = zeros(length(f_try), 3);  velocitiesp_try = coordsp_try;
        
    % helper functions
    r     = @(th) k0*exp(k1*sin(k2*th+phi));    
    thdot = @(th) sqrt( muC./r(th).^3 .* ...
        (1 + k1^2*k2^2*cos(k2*th+phi).^2 + k1*k2^2*sin(k2*th+phi)) );
    rdot  = @(th) thdot(th).*k1*k2*cos(k2*th+phi).*r(th);
    
    % initial coordinates & velocities 
    [x, y] = pol2cart(0, r(0));          [xdot, ydot] = pol2cart(thdot(0), rdot(0));
    coordsp_try(1, :) = [x,y,0]*R;       velocitiesp_try(1, :) = [xdot,ydot,0]*R;
    
    % kill warnings caused by QUADGK()
    warnstate_quadgk = warning('off', 'MATLAB:quadgk:NonFiniteValue');
    
    % loop through all elements (every theta-step needs an integration)
    for i = 1:numel(f_try)-1
        % rename angles for clarity
        th0 = f_try(i);  thend = f_try(i+1);        
        % get time (with Gauss/Kronrod quadrature)
        t_try(i+1) = t_try(i) + quadgk(@tff, th0, thend)/sqrt(muC);        
        % compute coordinates & velocities at this new [theta]
        [x, y] = pol2cart(thend, r(thend));  [xdot, ydot] = pol2cart(thdot(thend), rdot(thend)); 
        coordsp_try(i+1, :) = [x,y,0]*R;     velocitiesp_try(i+1, :) = [xdot,ydot,0]*R;  
    end
    
    % reset warnings to original setting
    warning(warnstate_quadgk); 
    
    % make cubic interpolations for every element in the statevector
    % (Every MATLAB version has SPLINE(), so it's safe to use it here)
    xcoefs = spline(t_try, coordsp_try(:, 1));   xdotcoefs = spline(t_try, velocitiesp_try(:, 1)); 
    ycoefs = spline(t_try, coordsp_try(:, 2));   ydotcoefs = spline(t_try, velocitiesp_try(:, 2));    
    zcoefs = spline(t_try, coordsp_try(:, 3));   zdotcoefs = spline(t_try, velocitiesp_try(:, 3));
    
    % In order to have realistic computation times for 100.000 MP's (or more), we
    % MUST convert this to constant-timestep cubic splines (use 1 day step):
    times = 0:86400:timespan;
    
    % make new cubic splines objects
    xcoefs = spline(times, ppval(xcoefs,times));  xdotcoefs = spline(times, ppval(xdotcoefs,times)); 
    ycoefs = spline(times, ppval(ycoefs,times));  ydotcoefs = spline(times, ppval(ydotcoefs,times)); 
    zcoefs = spline(times, ppval(zcoefs,times));  zdotcoefs = spline(times, ppval(zdotcoefs,times));  
    
    % save breakpoints 
    breaks = xcoefs.breaks;
    
    % arrange coefficients so that extraction from memory can be done at
    % peak efficiency (page-wise once, row-wise 6 times)
    spline_coefs = [reshape(xcoefs.coefs.'   , 4,1,[]), reshape(ycoefs.coefs.'   , 4,1,[]),...
                    reshape(zcoefs.coefs.'   , 4,1,[]), reshape(xdotcoefs.coefs.', 4,1,[]),...
                    reshape(ydotcoefs.coefs.', 4,1,[]), reshape(zdotcoefs.coefs.', 4,1,[])];
                
    %% INITIALIZE SECONDARIES
    
    % extract values from the secondaries
    as = secs(:, 1);   Omegas = secs(:, 4);
    es = secs(:, 2);   omegas = secs(:, 5);
    Is = secs(:, 3);   f0s    = secs(:, 6);
    
    % initialize some more
    warning_state = warning('off', 'MATLAB:divideByZero'); % just annoying!
    timespan = tend - t0;                   % total time window
    ns       = sqrt(muC./abs(as).^3);       % mean motions
    Ts       = 2*pi./ns;                    % periods
    M0s      = mod(etheta2M(es, f0s), 2*pi);% initial mean anomalies 
    Mends    = M0s + ns*timespan;           % final mean anomalies     
    fends    = eM2theta(es, Mends);         % final true anomalies 
    
    % initial & final positions & velocities 
    x0s   = kep2cart([as,es,Is,Omegas,omegas,f0s]  , muC, 'theta'); 
    xends = kep2cart([as,es,Is,Omegas,omegas,fends], muC, 'theta'); 
    
    % compute hte minimum and maximum distances to the central body
    [perip, app] = minmax_distances_exposin;
    [peris, aps] = minmax_distances(x0s, xends, f0s, fends, as, es);
    
    % optional parameters
    if (nargin == 6) || (isempty(params)), params = []; end   
    % optional: output function
    if isfield(params, 'outputFcn') && isa(params.outputFcn, 'function_handle')
        have_outputFcn = true;
        outputFcn = params.outputFcn;
    else
        have_outputFcn = false;
    end
  
    % optional: constant or distance-dependent threshold
    if isfield(params, 'threshold')
        % threshold depends on apocenter of primary
        if length(params.threshold) == 3
            threshold = params.threshold(1) + params.threshold(2)*app/params.threshold(3);
        % use constant threshold 
        else
            threshold = params.threshold; 
        end
    end  
    
    % evaluate output function
    if have_outputFcn
        % collect values
        state = 'initialization';
        values.candidates  = size(es,1);
        values.rejected    = 0;
        values.accepted    = 0;
        values.coplanar    = 0;
        values.noncoplanar = 0;
        % evaluate function
        outputFcn(state, values);
    end
    
    %% APSES/MIN-MAX-DISTANCE CHECK
    
    % if the pericenter (or minimum distance) of the primary lies above the apocenter 
    % of the secondary, or if the apocenter of the primary (or the maximum distance) 
    % lies below the pericenter of the secondary, the two can never encounter
    good_inds = all_inds((perip - aps <= threshold) & (peris - app <= threshold));
    
    % evaluate output function
    if have_outputFcn
        % collect values
        state = 'apsescheck';
        values.candidates  = size(es,1);
        values.rejected    = size(es,1) - numel(good_inds);
        values.accepted    = numel(good_inds);
        values.coplanar    = 0;
        values.noncoplanar = 0;
        % evaluate function
        outputFcn(state, values);
    end
    
    % they ALL might be unreachable. Exit in that case
    if isempty(good_inds), warning(warning_state); return, end

    %% MAIN CHECK
    
    % prune initial values of the secondaries
    x0s = x0s(good_inds, :);     
    Ts  = Ts(good_inds);
    
    % number of MP's to try
    num_canditates = size(x0s,1);
    
    % really the only way to determine the minimum distance from an MP to
    % an ExopSin, is by interpolating a high-order Chebychev polynomial
    % through a number of intermediate distances. Some experimentation has
    % shown that 10 such evaluations suffices to get a very rough initial
    % estimate while keeping the method reasonably fast, and 5 more 
    % evaluations provides an acceptable accuracy
    for i = 1:num_canditates
        
        % evaluate output function
        % NOTE: don't evaluate too often, it'll produce quite some overhead
        % NOTE: Use (ii-1) to make sure the first iteration is evaluated      
        if have_outputFcn && (mod(i-1, 5000) == 0)
            % collect values
            state = 'overlapcheck';
            values.candidates  = num_canditates;
            values.accepted    = nnz(reachable);
            values.rejected    = i - values.accepted;
            values.iteration   = i;
            stop = outputFcn(state, values);
            % when stop evaluates to [true], exit
            if stop, return, end
        end
        
        % find roots of the time derivative of the distance, using a
        % 10-node Chebychev interpolation
        t_roots = FindRealRoots(@(x) Rrrel(x,i), t0, tend, 10, true);
        % get the associated distances, and find the minimum        
        [dummy, distances] = Rrrel(t_roots,i);%#ok
        % find the root that yields the minimum distance
        [minimum_distance, which_t] = min(distances);
        
        % evaluate the function 5 more times to find a more
        % accurate approximation
        % NOTE: the inaccuracy introduced by using only 10 nodes might 
        % give rise to non-existent roots. To prevent these cases, we 
        % have to loop through all the minimal distances until a valid 
        % root was found
        while ~isempty(distances)
            better_root = FindRealRoots(@(x)Rrrel(x,i), ...
                max(t_roots(which_t)-Ts(i)/5, t0), ...
                min(t_roots(which_t)+Ts(i)/5, tend), 5, true);
            if isempty(better_root)  
                t_roots   = t_roots(distances ~= minimum_distance);
                distances = distances(distances ~= minimum_distance);                
                [minimum_distance, which_t] = min(distances);
            else break
            end
        end 
        % they ALL might have failed; then, just give up on this MP 
        if isempty(distances), continue, end
        % find the root that yields the minimum distance
        [dummy, minimum_distance, speed] = Rrrel(better_root,i);%#ok
        % there might *still* be multiple roots
        [minimum_distance, which_t] = min(minimum_distance);
        better_root = better_root(which_t); 
        
        % save results  
        if (minimum_distance <= threshold)
            min_dist       (good_inds(i)) = minimum_distance;
            rel_speed      (good_inds(i)) = speed(which_t);
            encounter_times(good_inds(i)) = better_root;
            reachable      (good_inds(i)) = true;
        end
    end % main loop  
    
    % evaluate output function
    % (this makes sure also the LAST iteration is evaluated)
    if have_outputFcn
        % collect values
        state = 'overlapcheck';
        values.candidates  = num_canditates;
        values.accepted    = nnz(reachable);
        values.rejected    = num_canditates - values.accepted;
        values.iteration   = size(es,1);
        % evaluate function
        outputFcn(state, values);
    end
    
    % scale encounter times back to days
    encounter_times = encounter_times/86400;
    
    % reset warnings
    warning(warning_state); 

    %% HELPER FUNCTIONS
     
    % compute instantaneous distance, and time-derivative thereof
    function [R, distance, speed, statep,states] = Rrrel(t, index) 
        % compute new states of the secondaries 
        statep = zeros(numel(t), 6);
        % (Code copy-pasted from  PPVAL() and adjusted for speed)
        for ii = 1:numel(t)
            % compute breakpoint (the cubic splines interpolation was 
            % forced to have constant time steps for THIS reason: 
            breakpoint = max(1, floor((t(ii)-t0)/86400));
            % shift times to proper evaluation sites (all breakpoints are the same)
            t2 = t(ii) - t0 - breaks(breakpoint);
            % calculate coordinates and speeds
            statep(ii, :) = [t2*t2*t2, t2*t2, t2, 1] * spline_coefs(:,:,breakpoint);
        end        
        % compute new states of the secondaries        
        states = progress_orbit(t(:)-t0,x0s(index,:), muC, 'seconds');
        % rename
        rrel = statep(:,1:3) - states(:, 1:3);
        Vrel = statep(:,4:6) - states(:, 4:6);
        % time derivative
        R = 2*sum(rrel.*Vrel,2);
        % distance
        distance = sqrt(sum(rrel.*rrel,2));
        % relative speed
        speed = sqrt(sum(Vrel.*Vrel,2));
    end % Rrrel
    
    % time of flight integrand
    function tof = tff(th)              
        % some substitutions
        k2thpphi = k2*th + phi;            
        s        = sin(k2thpphi);           c     = cos(k2thpphi);
        k12s     = k1*k2*k2*s;              tangm = k1*k2*c;
        tanfc    = tangm.^2 + k12s + 1;     r32   = k0^(3/2)*exp(3*k1*s/2);        
        % equation for transfer time
        tof = r32.*sqrt(abs(tanfc)); % abs is just to be sure        
    end % tof
    
    % Compute minimum and maximum distance to the central body for the ExpoSin
    function [minimum_distance, maximum_distance] = minmax_distances_exposin
                
        % compute maximum and minimum values of
        % sin(k2*th + phi), with 0 < th < dth
        if (dth >= 2*pi/k2) || (mod(phi,2*pi) < pi/2 && mod(k2*dth + phi, 2*pi) > pi/2)
            maxsin = +1;
        else
            maxsin = max(sin(phi), sin(k2*dth + phi));
        end
        if (dth >= 2*pi/k2) || (mod(phi,2*pi) < 3*pi/2 && mod(k2*dth + phi, 2*pi) > 3*pi/2)
            minsin = -1;
        else
            minsin = min(sin(phi), sin(k2*dth + phi));
        end
        
        % determine minimum and maximum values of
        % k1 * sin(k2*th + phi), with 0 < th < dth
        if (k1 > 0)
            mink1sin = k1*minsin;   maxk1sin = k1*maxsin;
        else
            mink1sin = k1*maxsin;   maxk1sin = k1*minsin;
        end
        
        % minimum and maximum distances
        minimum_distance = k0*exp(mink1sin);
        maximum_distance = k0*exp(maxk1sin);
        
    end % min/max distances

    % compute minimum and maximum distances to the central body for the secondaries
    function [minimum_distance, maximum_distance]  = minmax_distances(...
            r1vec, r2vec, f0, fend, a, e)
        
        %% initialize
                
        % compute r1, r2
        r1 = sqrt(sum(r1vec.*r1vec,2));
        r2 = sqrt(sum(r2vec.*r2vec,2));
        
        % default - minimum of r1,r2
        minimum_distance = min([r1,r2],[],2);
        maximum_distance = max([r1,r2],[],2);
        
        % compute dth 
        dth = fend - f0;
        
        % was the longway used or not?
        longway = abs(dth) > pi;
        
        % NOW put both [theta]'s in [-pi, pi]
        f0   = mod(f0, 2*pi);     f0(f0 > pi)     = f0(f0 > pi) - 2*pi;
        fend = mod(fend, 2*pi);   fend(fend > pi) = fend(fend > pi) - 2*pi;
        
        % calculate apses
        pericenter = a.*(1-e);
        apocenter  = zeros(size(pericenter));
        apocenter(e >= 1) = inf; % parabolic/hyperbolic case
        apocenter(e <  1) = a(e < 1).*(1+e(e < 1)); % elliptic case
        
        %% obvious cases 
        
        % both apses are obviously traversed
        minimum_distance(dth > 2*pi) = pericenter(dth > 2*pi);
        maximum_distance(dth > 2*pi) = apocenter(dth > 2*pi);
        
        %% less obvious cases
        
        % points 1&2 are on opposite sides of the symmetry axis -- minimum
        % and maximum distance depends both on the value of [dth], and both
        % [theta1] and [theta2]
        opposites = (f0.*fend < 0);
        
        % if |th1| + |th2| = turnangle, we know that the pericenter was passed
        change_mindist = opposites & (abs(f0)+abs(fend) == dth);
        minimum_distance(change_mindist) = pericenter(change_mindist);
        
        % this condition can only be false for elliptic cases, and when it is 
        % indeed false, we know that the orbit passed apocenter
        change_maxdist = opposites & ~change_mindist;
        maximum_distance(change_maxdist) = apocenter(change_maxdist);
        
        % points 1&2 are on the same side of the symmetry axis. Only if the
        % long-way was used are the min. and max. distances different from
        % the min. and max. values of the radii (namely, equal to the apses)            
        change_sames = ~opposites & longway;        
        minimum_distance(change_sames) = pericenter(change_sames);
        maximum_distance(change_sames & (e < 1)) = apocenter(change_sames & (e < 1));
        
    end
    
end % minimum distance exposins
