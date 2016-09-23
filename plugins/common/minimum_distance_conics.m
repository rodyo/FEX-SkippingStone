function [reachable, encounter_times, min_dist, rel_speed] = minimum_distance_conics(...
        pri, secs, t0, tend, muC, params)
% MINIMUM_DISTANCE          Calculates the minimum distance between an
%                           arbitrary conic section, and an arbitrary
%                           number of ellipses
%
% USAGE:
% [reachable, encounter_times, min_dist] = ...
%          MINIMUM_DISTANCE(pri, secs, threshold, t0, tend, muC)
%
% INPUT ARGUMENTS:
%     pri : primary conic section to which the distances are 
%           to be computed (usually the S/C trajectory). 
%     sec : all secondary conic sections.
%
%           Both can be ellipses, parabolae or hyperbolae, to 
%           be given in the format
%
%              pri = [a, e, I, Omega, omega, f0]
%
%           and
%
%              sec = [a(1), e(1), I(1), Omega(1), omega(1), f0(1); 
%                     a(2), e(2), I(2), Omega(2), omega(2), f0(2);
%                     ...
%                     a(N), e(N), I(N), Omega(N), omega(N), f0(N)];
%
%           where a is the semi-major axis, e the eccentricity, I the
%           inclination, Omega the longitude of the ascending node,
%           omega the argument of pericenter and f0 the true anomaly
%           at t = t0.
%
% t0/tend : start/end time [days]. 
%     muC : standard gravitational parameter of the central body [km3s-2].
%  params : parameter-value pairs. Valid combinations are:
%
%           'outputFcn', function_handle    
%              define an output function. defaults to [] (none).
%
%           'threshold', value(s)           
%              distance below which the corresponding secondary is
%              considered to be 'close'. Can be a single value [km], 
%              or a 3-element vector of values, in which case a
%              variable theshold [T] will be used based on the
%              periapsis of the corresponding secondary, via the
%              formula:
%
%                   threshold = T(1) + T(2)*a(sec)/T(3)
%
%              with a(sec) the pericenter distance of the secondary. 
%
% OUTPUT ARGUMENTS:
%
%      reachable : logical indices to [sec], indicating which 
%                  secondaries come closer than the given threshold 
%                  (within [t0 tend]).
%encounter_times : Corresponding encounter times
%       min_dist : corresponding minimum distances [km]
%      rel_speed : relative speed at the closest encounter [km/s]
%
%
% NOTE: Generally this is a very CPU intensive computation. The code has
% been heavily optimized for speed, but remains coded in M. Ideally, this
% whole function is either translated into an EML compatible format, or
% re-coded in C/C++ to improve execution speed. At the time of writing 
% I didn't have the energy left to do this (as this is quite a bit of 
% work), so this is left for some other time. 
%                                   - Rody Oldenhuis
%

% This is an implementation of the methods described in 
%
%  Hoots F.R., Crawford L.L., and Roehrich R.L. "An analytic method to 
%  determine future close approaches between satellites". Celestial Mechanics, 
%  33:143–158, 1984. DOI 0008-8714/84/0332.
%
% with some modifications and optimizations as described in 
%
%  "Trajectory Optimization of a Mission to the Solar Bow Shock and Minor
%  Planets", Chapter 11.4 and onward. MSc. thesis by Rody P.S. Oldenhuis, 
%  Delft University of Technology, Workgroup "Astrodynamics and Satellite 
%  Missions", Feb/2010.

% Author: Rody P.S. Oldenhuis
% M.Sc. student in the Delft University of Technology
%
% Partly performed at ISAS/JAXA, 
% October 2008 - April 2009
% Last updated for thesis research, 31/Oct/2009

      
    %% initialize
    
    % parameters & constants   
    small_e   = 0.3; % define what a "small" eccentricity is
    threshold = 0.01 * 1.495978706910000e+008; % default threshold; 0.01 AU        
    warning_state = warning('off', 'MATLAB:divideByZero'); % just annoying!
    
    % errortrap (ALL arguments must be given!)
    error(nargchk(5,6,nargin));
    
    % times are given in days; convert to seconds
    % (everything else is given in seconds)
    t0   = t0 * 86400;
    tend = tend * 86400;
    
    % extract values from the primary
    ap = pri(1);       Omegap = pri(4);
    ep = pri(2);       omegap = pri(5);
    Ip = pri(3);       f0p    = pri(6);
        
    % extract values from the secondaries
    as = secs(:, 1);   Omegas = secs(:, 4);
    es = secs(:, 2);   omegas = secs(:, 5);
    Is = secs(:, 3);   f0s    = secs(:, 6);

    % compute unit vectors wp and ws (Eq. 3)
    % NOTE: unlike described in the paper, these vectors are not ACTUALLY
    % vectors perpendicular to the orbital plane; these vectors just
    % produce the correct value for sin(IR) when taking the cross-product
    wp = [sin(Omegap) *sin(Ip), cos(Omegap) *sin(Ip), cos(Ip)]; % 1x3 vector
    ws = [sin(Omegas).*sin(Is), cos(Omegas).*sin(Is), cos(Is)]; % Nx3 matrix
    
    % compute IRs (Eq. 4)
    wp  = wp(ones(size(ws,1),1), :); 
    K1  = [ws(:,2).*wp(:,3) - ws(:,3).*wp(:,2)
           ws(:,3).*wp(:,1) - ws(:,1).*wp(:,3)
           ws(:,1).*wp(:,2) - ws(:,2).*wp(:,1)];
    IRs = asin(sqrt(sum(K1.*K1, 2)));
    
    % clear unneeded variables
    clear w0 wp ws K1;

    % compute Deltapri (Eqs. 8 & 9)
    cosDeltap = (sin(Ip)*cos(Is) - cos(Ip)*sin(Is).*cos(Omegap-Omegas));
    sinDeltap = (sin(Is).*sin(Omegap-Omegas));
    Deltap    = atan2(sinDeltap, cosDeltap);
    
    % clear unneeded variables
    clear cosDeltap sinDeltap

    % Compute Deltasec (Eqs. 10 & 11)
    cosDeltas = (sin(Ip)*cos(Is).*cos(Omegap-Omegas) - sin(Is)*cos(Ip));
    sinDeltas = (sin(Ip)*sin(Omegap-Omegas));
    Deltas    = atan2(sinDeltas, cosDeltas);
    
    % clear unneeded variables
    clear cosDeltas sinDeltas

    % substitution constants (Eq. 17)
    axp = ep*cos(omegap - Deltap);      axs = es.*cos(omegas - Deltas);
    ayp = ep*sin(omegap - Deltap);      ays = es.*sin(omegas - Deltas);

    % mean motions & periods
    np = sqrt(muC /abs(ap)^3);     Tp = 2*pi /np *(ep < 1);  % zero when hyperbolic
    ns = sqrt(muC./abs(as).^3);    Ts = 2*pi./ns.*(es < 1);  % zeor when hyperbolic 
    
    % total time window
    timespan = tend - t0;
    
    % initial & final mean anomalies     
    if (ep < 1)
        M0p = mod(etheta2M(ep, f0p), 2*pi);
    else
        M0p = etheta2M(ep, f0p);
    end,                                  Mendp = M0p + np*timespan;
    M0s = mod(etheta2M(es, f0s), 2*pi);   Mends = M0s + ns*timespan;
    
    % convert final [M] to [theta]
    fendp = eM2theta(ep, Mendp);
    fends = eM2theta(es, Mends);
    
    % initial & final positions & velocities    
    x0p   = kep2cart([ap,ep,Ip,Omegap,omegap,f0p],   muC, 'theta'); % NOTE: using [theta]
    x0s   = kep2cart([as,es,Is,Omegas,omegas,f0s],   muC, 'theta'); % is much faster than 
    xendp = kep2cart([ap,ep,Ip,Omegap,omegap,fendp], muC, 'theta'); % using 'M'
    xends = kep2cart([as,es,Is,Omegas,omegas,fends], muC, 'theta'); %
    
    % clear unneeded variables
    clear Omegas Omegap Is Ip;
    
    % compute the minimum and maximum distances to the central body
    [perip, app] = minmax_distances(x0p, xendp, f0p, fendp, ap, ep);
    [peris, aps] = minmax_distances(x0s, xends, f0s, fends, as, es);
    
    % clear unneeded variables
    clear f0p f0s fendp fends;
                
    % initial and final distances and relative speeds
    dist0    = sqrt(sum((bsxfun(@minus,   x0p(:,1:3),   x0s(:,1:3) )).^2, 2));
    distend  = sqrt(sum((bsxfun(@minus, xendp(:,1:3), xends(:,1:3) )).^2, 2));
    speed0   = sqrt(sum((bsxfun(@minus,   x0p(:,4:6),   x0s(:,4:6) )).^2, 2));
    speedend = sqrt(sum((bsxfun(@minus, xendp(:,4:6), xends(:,4:6) )).^2, 2));
        
    % clear unneeded variables
    clear xendp xends;
    
    % how many complete periods fit in the time window
    complete_Periods_pri = ceil((Mendp - M0p)/2/pi) *(ep < 1); % zero when hyperbolic
    complete_Periods_sec = ceil((Mends - M0s)/2/pi).*(es < 1);
    
    % initialize output values (it's easiest to be pessimistic)
    reachable       = false(size(secs,1),1);
    encounter_times = NaN(size(secs,1),1);
    min_dist        = inf(size(secs,1),1);
    rel_speed       = min_dist;
    
    % optional: minimum/maximum distance of primary to the central body
    if (nargin == 5) || (isempty(params)), params = []; end   
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
            threshold = params.threshold(1) + ...
                params.threshold(2)*min(app,30*150e6)/params.threshold(3);
        % use constant threshold 
        else
            threshold = params.threshold; 
        end
    end  
    
    % clear unneeded variables
    clear('params');
        
    % compute square of the threshold only once
    threshold2 = threshold^2;  
    
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
        
    %% Apses/min-max-distance check
    
    % if the pericenter (or minimum distance) of the primary lies above the apocenter 
    % of the secondary, or if the apocenter of the primary (or the maximum distance) 
    % lies below the pericenter of the secondary, the two can never encounter
    good_inds = find( (perip-aps <= threshold) & (peris-app <= threshold) );
    
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
    
    %% Time prefilter & compute minimum distances
    
    % For those still in the game, prune with the time prefilter, and
    % compute the initial values for the time-dependent min. distance
            
    % prune values
    axp    = axp(good_inds);        axs    = axs(good_inds);
    ayp    = ayp(good_inds);        ays    = ays(good_inds);
    as     = as(good_inds);         es     = es(good_inds);        
    IRs    = IRs(good_inds);        ns     = ns(good_inds);
    Deltap = Deltap(good_inds);     Deltas = Deltas(good_inds);        
    omegas = omegas(good_inds);     M0s    = M0s(good_inds);    
    x0s    = x0s(good_inds, :);     Ts     = Ts(good_inds);    
    aps    = aps(good_inds);        peris  = peris(good_inds);
    complete_Periods_sec = complete_Periods_sec(good_inds);
        
    % compute Qp & Qs (Eq. 28)
    alphap = ap *(1 - ep ^2).*sin(IRs);      
    alphas = as.*(1 - es.^2).*sin(IRs);
    Qp = alphap.*(alphap - 2*threshold*ayp) - (1 - ep ^2)*threshold2;
    Qs = alphas.*(alphas - 2*threshold*ays) - (1 - es.^2)*threshold2;

    % clear unneeded variables
    clear IRs;
    
    % compute both cosines of uRp and uRs (Eq. 27)
    cosuRp1 = (-threshold2*axp + (alphap - threshold*ayp).*sqrt(Qp)) ./ ...
              (alphap.*(alphap - 2*threshold*ayp) + threshold2*ep ^2);
    cosuRp2 = (-threshold2*axp - (alphap - threshold*ayp).*sqrt(Qp)) ./ ...
              (alphap.*(alphap - 2*threshold*ayp) + threshold2*ep ^2);
    cosuRs1 = (-threshold2*axs + (alphas - threshold*ays).*sqrt(Qs)) ./ ...
              (alphas.*(alphas - 2*threshold*ays) + threshold2*es.^2);
    cosuRs2 = (-threshold2*axs - (alphas - threshold*ays).*sqrt(Qs)) ./ ...
              (alphas.*(alphas - 2*threshold*ays) + threshold2*es.^2);  
          
    % clear unneeded variables     
    clear axp axs ayp ays
          
    % some of these must be computed with the coplanar method. 
    % Find which ones they are
    coplanar_inds = (Qp < 0) | abs(cosuRp1) > 1 | abs(cosuRp2) > 1 | ...
                    (Qs < 0) | abs(cosuRs1) > 1 | abs(cosuRs2) > 1;     
                
    % clear unneeded variables     
    clear Qp Qs;
                
    % rename indices (it's just more readible)
    nc = ~coplanar_inds;                   % indices to non-coplanar cases
    cp =  coplanar_inds & (es <= small_e); % indices to low-eccentricity coplanar cases
    
    % clear the old ones
    clear coplanar_inds;
    
    % number of MP's per case (for in the output function)
    num_coplanar    = nnz(cp);    num_canditates = size(es,1);
    num_noncoplanar = nnz(nc);
    
    % initialize some variables
    uRp     = zeros(size(es,1),4);   uRs     = uRp;  
    good_tp = uRp;                   good_ts = uRp;
    
    % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    % initialize non-coplanar cases
    % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if any(nc)
        
        % compute angular windows (Eq. 29) 
        acosuRp1   = acos(cosuRp1(nc));   acosuRs1   = acos(cosuRs1(nc));
        acosuRp2   = acos(cosuRp2(nc));   acosuRs2   = acos(cosuRs2(nc));
        uRp(nc, 1) = +acosuRp1;           uRs(nc, 1) = +acosuRs1;
        uRp(nc, 2) = -acosuRp1;           uRs(nc, 2) = -acosuRs1;
        uRp(nc, 3) = +acosuRp2;           uRs(nc, 3) = +acosuRs2;
        uRp(nc, 4) = -acosuRp2;           uRs(nc, 4) = -acosuRs2;

        % compute corresponding true anomalies (Eq. 14)
        if (ep < 1)
            fps = mod(bsxfun(@plus, uRp(nc,:), Deltap(nc) - omegap), 2*pi);
        else
            fps = bsxfun(@plus, uRp(nc,:), Deltap(nc) - omegap);
            if any(fps(:) > pi)
                fpslpi = find(fps > pi);
                fps(fpslpi) = fps(fpslpi) - 2*pi;
            end
        end
        fss = mod(bsxfun(@plus, uRs(nc,:), Deltas(nc) - omegas(nc)), 2*pi); 

        % convert to times via Kepler's equation         
        new_Mp = etheta2M(repmat(ep, size(fps,1),4), fps) - M0p;        
        new_Ms = etheta2M(repmat(es(nc), 1,4), fss) - repmat(M0s(nc), 1,4);        
        good_tp(nc,:) = new_Mp/np + t0;
        good_ts(nc,:) = bsxfun(@rdivide,new_Ms,ns(nc)) + t0;
               
        % clear unneeded variables
        clear acosuRp1 acosuRs1 acosuRp2 acosuRs2 fss fps;
        clear cosuRp1 cosuRs1 cosuRp2 cosuRs2 new_Mp new_Ms;
        
    end
    
    % clear unneeded variables
    clear uRp uRs Deltap Deltas ns omegas;
    
    % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    % initialize low-eccentricity coplanar cases
    % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    if any(cp)
        
        % use the averages of the apses as the circular radius
        rss = (aps(cp) + peris(cp))/2;

        % compute true and mean anomalies at which the primary will 
        % intersect the secondary         
        theta01 = +real(acos((ap*(1-ep^2)./rss - 1)/ep));   M01 = etheta2M(ep, theta01);
        theta02 = -real(acos((ap*(1-ep^2)./rss - 1)/ep));   M02 = etheta2M(ep, theta02);
        
        % the times of these intersections are now easy to determine:
        t0p1 = (M01 - M0p)/np + t0;
        t0p2 = (M02 - M0p)/np + t0;    

        % convert these initial times to a time window. These time windows are
        % centered around these initial times, and are 1/5 of the period of the
        % secondary wide
        good_tp(cp,:) = sort([t0p1-Ts(cp)/5, t0p1+Ts(cp)/5, t0p2-Ts(cp)/5, t0p2+Ts(cp)/5], 2);
        good_ts(cp,:) = good_tp(cp,:);
        
        % clear unneeded variables
        clear t0p1 t0p2 M01 M02 theta01 theta02 np rss;
        
    end
    
     % clear unneeded variables
    clear M0s M0p aps peris;
    
    % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    % start solution method
    % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     
    % initialize
    MPs           = size(es,1);
    can_reach     = false(MPs,1);
    mindistance   = inf(MPs,1);
    relativespeed = mindistance;
    enctimes      = NaN(MPs,1);
    
    % nested loop
    for ii = 1:MPs
        
        % initialize
        mindist       = inf;
        relspeed      = inf;
        encountertime = NaN;
        
        % evaluate output function
        % NOTE: don't evaluate too often, it'll produce quite some overhead
        % NOTE: Use (ii-1) to make sure the first iteration is evaluated
        if have_outputFcn && (mod(ii-1, 5000) == 0)
            % collect values
            state = 'overlapcheck';
            values.candidates  = num_canditates;
            values.accepted    = nnz(reachable);
            values.rejected    = ii - values.accepted;
            values.coplanar    = num_coplanar;
            values.noncoplanar = num_noncoplanar;
            values.iteration   = ii;
            stop = outputFcn(state, values);
            % when stop evaluates to [true], exit
            if stop, return, end
        end
               
        % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        % Non-coplanar or low-eccentricity coplanar case: Check for overlap
        % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        %
        % There are two possible cases:
        %
        %  t01--------te1                         t01---te1
        %        t02-------------te2   or     t02-----------te2
        %
        % The condition for overlap for both of these cases is simply
        %
        %                 (t02 <= te1) & (te2 >= t01)
        %
        % The condition for overlap has to be computed by comparing each
        % value in the time vector of the secondary, to each element in the
        % time vector of the primary.
        %
        % Instead of using Newton's method, it is more efficient (and
        % robust) to use a Chebychev-interpolation on the whole interval
        % and find all of its roots. Doing it that way there is no chance
        % the routine gets stuck on a maximum or in an infinite loop, and
        % the number of function evaluations is always constant.
        if nc(ii) || cp(ii)
            
            % add complete multiples to this primary
            % NOTE: avoid repmat to reduce overhead
            add_this0 = 0:min(100,complete_Periods_pri); % EML requires bounding the variable
            add_this1  = add_this0(ones(4,1), :) * Tp;
            all_good_tps0 = good_tp(ii,:).';
            all_good_tps1 = all_good_tps0(:,...
                ones(1,max(1, min(100,complete_Periods_pri+1))));
            all_good_tps  = sort(all_good_tps1(:).' + add_this1(:).');
            
            % extract windows
            t01 = all_good_tps(1:2:end); % start of each window
            te1 = all_good_tps(2:2:end); % end of each window
            
            % handle constraints
            too_early  = (t01 < t0  ) & (te1 < t0  );  t01 = t01(~too_early);  te1 = te1(~too_early);
            too_late   = (t01 > tend) & (te1 > tend);  t01 = t01(~too_late);   te1 = te1(~too_late);
            half_early = (t01 < t0  ) & (te1 > t0  );  t01(half_early) = t0;
            half_late  = (t01 < tend) & (te1 > tend);  te1(half_late)  = tend;
            
            % number of windows for the primary
            num_window_pris = numel(t01);
            
            % add complete multiples of the period to this secondary
            % NOTE: avoid repmat to reduce overhead
            add_this2 = 0:min(100,complete_Periods_sec(ii));
            add_this3 = add_this2(ones(4,1), :) * Ts(ii);
            all_good_tss2 = good_ts(ii, :).';
            all_good_tss3 = all_good_tss2(:, ones(1,max(1, min(100,complete_Periods_sec(ii)+1))));
            all_good_tss = sort(all_good_tss3(:).' + add_this3(:).');
            
            % extract windows
            t02 = all_good_tss(1:2:end);  % start of each window
            te2 = all_good_tss(2:2:end);  % end of each window
            
            % handle constraints
            too_early  = (t02 < t0  ) & (te2 < t0  );  t02 = t02(~too_early);   te2 = te2(~too_early);
            too_late   = (t02 > tend) & (te2 > tend);  t02 = t02(~too_late);    te2 = te2(~too_late);
            half_early = (t02 < t0  ) & (te2 > t0  );  t02(half_early) = t0;
            half_late  = (t02 < tend) & (te2 > tend);  te2(half_late)  = tend;
            
            % number of windows for the secondary
            num_window_secs = numel(t02);
            
            % they might fall entirely outside of the time span;
            % continue in that case
            if isempty(t01) || isempty(te1) || isempty(t02) || isempty(te2), continue, end
            
            % compute minima and maxima only once
            maxte1 = max(te1);      maxte2 = max(te2);
            mint01 = min(t01);      mint02 = min(t02);
            
            % only enter loops when needed
            if (maxte1 > mint02) && (maxte2 > mint01)
                % loop through each window of the primary
                for jj = 1:num_window_pris
                    % reduce number of iterations
                    if (t01(jj) > maxte2) || (te1(jj) < mint02), continue, end
                    % loop through each window of the secondaries
                    for kk = 1:num_window_secs
                        % reduce number of iterations
                        if (t02(kk) > maxte1) || (te2(kk) < mint01), break, end
                        % if true, compute the minimum distance and encounter time
                        if (t02(kk) <= te1(jj)) && (te2(kk) >= t01(jj))
                            % generate interval
                            a = max(t01(jj),t02(kk));   b = min(te1(jj),te2(kk));
                            % get distance in the middle of the interval. If it is
                            % larger than 25*threshold, the bodies will most
                            % likely be on opposite sides of their trajectory, and
                            % the actual root-finding can be skipped
                            %
                            % NOTE: This line will show up in the profiler
                            % as the most time-consuming line of all. But
                            % trust me, if you remove it, it'll be a LOT slower!
                            [~, try_distance] = Rrrel((a+b)/2, t0,x0p,x0s(ii,:),muC);
                            if (try_distance > 25*threshold), continue; end
                            % find the roots (if any) of the time-derivative of R
                            t = FindRrrelRoots(a, b, 5, t0,x0p,x0s(ii,:),muC);
                            % there might not be any roots; continue
                            if isempty(t), continue; end
                            % find the associated distances
                            [~,distances,speeds] = Rrrel(t, t0,x0p,x0s(ii,:),muC);
                            % there might also be multiple roots; find the one that
                            % yields the minimum distance
                            [minimum_distance, which_t] = min(distances);
                            % save its value if its less than the current minimum
                            if (minimum_distance < min_dist(good_inds(ii)))
                                mindist  = minimum_distance;
                                relspeed = speeds(which_t);
                                encountertime = t(which_t);
                            end
                        end % if overlap
                    end % time windows of secondaries (inner loop)
                end % time windows of primary (outer loop)
            end % if needed
            
        % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        % High-eccentricity coplanar case: Fit high-degree Chebychev polynomial
        % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=            
        else
            % evaluate the derivative at 25 nodes to find a descent
            % approximation to the root yielding the minimum
            t = FindRrrelRoots(t0, tend, 25, t0,x0p,x0s(ii,:),muC);
            % if there aren't any real roots, continue
            if isempty(t), continue,  end
            % retreive the associated distances
            [~, distances] = Rrrel(t, t0,x0p,x0s(ii,:),muC);
            % find the root that yields the minimum distance
            [minimum_distance, which_t] = min(distances);
            % evaluate the function 5 more times to find a more
            % accurate approximation
            % NOTE: the inaccuracy introduced by using only 25 nodes might
            % give rise to non-existent roots. To prevent these cases, we
            % have to loop through all the minimal distances until a valid
            % root was found
            better_root = FindRrrelRoots(...
                t(which_t)-Ts(ii)/5, t(which_t)+Ts(ii)/5, 5, t0,x0p,x0s(ii,:),muC);
            while ~isempty(distances)
                if isempty(better_root)
                    t         = t(distances ~= minimum_distance);
                    distances = distances(distances ~= minimum_distance);
                    [minimum_distance, which_t] = min(distances);
                    try % [t] might be empty
                        better_root = FindRrrelRoots(...
                            t(which_t)-Ts(ii)/5, t(which_t)+Ts(ii)/5, 5, t0,x0p,x0s(ii,:),muC);
                    catch, break; %#ok
                    end
                else break
                end
            end
            % they ALL might have failed; then, just give up on this MP
            if isempty(distances) || isempty(better_root), continue, end
            % find the root that yields the minimum distance
            [~, minimum_distance0, speed] = Rrrel(better_root, t0,x0p,x0s(ii,:),muC);
            % there might *still* be multiple roots
            [minimum_distance, which_t] = min(minimum_distance0);
            better_root = better_root(which_t);
            % rename for consistency
            t = better_root; which_t = 1;
            % save results
            mindist       = minimum_distance;
            relspeed      = speed(which_t);
            encountertime = t(which_t);
            
        end % select method based on eccentricity of the secondary
        
        % check whether the minimum distance found is actually less
        % than the distance at the boundaries of the interval
        if (dist0(good_inds(ii)) < mindist)
            mindist  = dist0(good_inds(ii));
            relspeed = speed0(good_inds(ii));
            encountertime = t0;
        end
        if (distend(good_inds(ii)) < mindist)
            mindist  = distend(good_inds(ii));
            relspeed = speedend(good_inds(ii));
            encountertime = tend;
        end
        
        % only save values if the minimum distance thus found is less
        % than the threshold distance
        if (mindist <= threshold)
            can_reach(ii)     = true;
            mindistance(ii)   = mindist;
            relativespeed(ii) = relspeed;
            enctimes(ii)      = encountertime; 
        end
        
    end % main for
    
    % assign data to correct indices
    reachable(good_inds)       = can_reach;
    encounter_times(good_inds) = enctimes;
    min_dist(good_inds)        = mindistance;
    rel_speed(good_inds)       = relativespeed;
    
    % scale encounter times back to days
    encounter_times = encounter_times/86400;
    
    % reset warnings
    warning(warning_state); 
        
end % main function

%% subfunctions

% compute minimum and maximum distances to the central body
function [minimum_distance, maximum_distance] = minmax_distances(...
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
    f0 = mod(f0, 2*pi);        f0(f0 > pi) = f0(f0 > pi) - 2*pi;    
    fend = mod(fend, 2*pi);    fend(fend > pi) = fend(fend > pi) - 2*pi;    
    
    % calculate apses
    pericenter = a.*(1-e);
    apocenter  = zeros(size(pericenter));
    apocenter(e >= 1) = inf; % parabolic/hyperbolic case
    apocenter(e < 1) = a(e < 1).*(1+e(e < 1)); % elliptic case
    
    
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
        
end % min/max distances

% computes the instantaneous distance, its time derivative and the 
% instantaneous relative speed
function [R, distance, speed] = Rrrel(t, t0,x0p,x0s,muC) 
    % compute new states with progress_orbit()
    statep = progress_orbit(t(:)-t0,x0p, muC, 'seconds');
    states = progress_orbit(t(:)-t0,x0s, muC, 'seconds');
    % rename
    rrel = statep(:,1:3)-states(:, 1:3);
    Vrel = statep(:,4:6)-states(:, 4:6);
    % time derivative
    R = 2*sum(rrel.*Vrel,2);
    % distance
    distance = sqrt(sum(rrel.*rrel,2));
    % relative speed
    speed = sqrt(sum(Vrel.*Vrel,2));
end % Rrrel

% computes roots of the Rrrel-equation using Chebyshev polynomials
function Roots = FindRrrelRoots(a, b, n, t0,x0p,x0s,muC)

    % parse input and initialize.    
    if n <= 2, n = 3; end  % Minimum [n] is 3:  
    
    % some convenient variables
    bma = (b-a)/2;  bpa = (b+a)/2;   Roots = [];
    
    % Obtain the Chebyshev coefficients for the function
    %
    % Based on the routine given in Numerical Recipes (3rd) section 5.8;
    % calculates the Chebyshev coefficients necessary to approximate some
    % function over the interval [a,b]
    
    % initialize 
    c = zeros(1,n);  k=(1:n)';  y = cos(pi*((1:n)-1/2)/n); 
    % evaluate Rrrel on Chebychev nodes 
    f = Rrrel((y*bma)+bpa, t0,x0p,x0s,muC);
    
    % compute the coefficients
    for j=1:n, c(j)=(f(:).'*(cos((pi*(j-1))*((k-0.5)/n))))*(2-(j==1))/n; end       
        
    % coefficients may be [NaN] if [inf]
    % ??? TODO - it is of course possible for c(n) to be zero...
    if any(~isfinite(c(:))) || (c(n) == 0), return; end
        
    % Define [A] as the Frobenius-Chebyshev companion matrix. This is based
    % on the form given by J.P. Boyd, Appl. Num. Math. 56 pp.1077-1091 (2006).
    one = ones(n-3,1);
    A = diag([one/2; 0],-1) + diag([1; one/2],+1);
    A(end, :) = -c(1:n-1)/2/c(n);
    A(end,end-1) = A(end,end-1) + 0.5;
    
    % Now we have the companion matrix, we can find its eigenvalues using the
    % MATLAB built-in function. We're only interested in the real elements of
    % the matrix:
    eigvals = eig(A);  realvals = eigvals(imag(eigvals)==0);
        
    % if there aren't any real roots, return
    if isempty(realvals), return; end
    
    % Of course these are the roots scaled to the canonical interval [-1,1]. We
    % need to map them back onto the interval [a, b]; we widen the interval just
    % a tiny bit to make sure that we don't miss any that are right on the 
    % boundaries.    
    rangevals = nonzeros(realvals(abs(realvals) <= 1+1e-5));

    % also sort the roots
    Roots = sort(rangevals*bma + bpa);
            
end % FindRrrelRoots