%#eml
function [V1, V2, extremal_distances, endmass, exposin, R, exitflag] = ...
        lambert_low_exposins(r1vec, r2vec, tf, k2, N, M0, Isp, muC)
% LAMBERT_LOW                    Low-Thrust Lambert-targeter
%
%   Usage:
%   [V1, V2, exitflag, exposin, R] = LAMBERTLOW(r1, r2, tf, k2, N, mu_central)
%
%   LAMBERT_LOW() determines the exponential sinuisoid between [r1] = [x1 y1 z1] 
%   and [r2] = [x2 y2 z2], that takes a transfer time [tf] in days. [V1] and [V2] 
%   are the velocity vectors at [r1] and [r2], and [exposin] is the exponential 
%   sinusoid in the form [k0,k1,k2,phi,tf,N, dth,gamma1,gamma_m,gamma_M] that 
%   satisfies the constraints. The output argument [R] is the Euler-rotation 
%   matrix used to properly orient the 2-dimensional solution into the original 
%   3-dimensional geometry. 
%
%   In practice it is quite hard to find combinations of [k2] and [N] that make it 
%   possible to satisfy all constraints. If for some reason any of the constraints 
%   can not be met, all results will be NaN. The exact reason for why the 
%   procedure failed is embedded in the output argument [exitflag]:
%
%       exitflag = +1        Successful exit. 
%       exitflag = -1        The given problem has no solution.
%       exitflag = -2        The procedure terminated successfully, but one of 
%                            the constraints could not be met. 
%       exitflag = -3        The problem has a solution and all constraints
%                            are met, but a numerical under- or overflow
%                            was encountered. 
%       exitflag = -4        The procedure failed to find a root of the
%                            TOF-integral. 
%
%   Selection of the short or long way may be accomplished by setting a
%   *negative* value for those transfer times for which the *long* way is desired.
%   For all transfer times that are positive, the short way will be used. Similarly,
%   if multiple solutions appear to exist, the solution that will be returned depends
%   on the sign of [N]. A *negative* [N] will result in the left-branch, while a
%   *positive* value will result in the right-branch.
%
%   LAMBERT_LOW uses the method developped by D. Izzo in his 2006 paper,
%   which was based on the exponential sinusoids introduced by Petropoulos in his
%   PhD. dissertation.

% Author
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% last edited 05/Oct/2009    

% ADJUSTED FOR EML-COMPILATION 05/Oct/2009

    % parameters
    tol1       = 1e+1;        % 10 second accuracy on time-of-flight integral
    tol2       = 1e-6;        % constraint equation accuracy to be met
    longway    = sign(tf);    % short or long way depends in sign of [tf]
    leftbranch = sign(N);     % left or right branch depends on sign of [N]
    tf         = 86400*tf;    % convert time to seconds    

    % initial output is pessimistic
    exposin = NaN(1, 10);    
    V1      = NaN(1, 3);
    V2      = V1;
    R       = NaN(3);
    endmass = NaN;
    extremal_distances = NaN(1,2);

    % initial values
    r1      = sqrt(r1vec*r1vec.');    % magnitude of r1
    r2      = sqrt(r2vec*r2vec.');    % magnitude of r2
    dotprod = (r1vec*r2vec.')/r1/r2;  % dot product of the two
    ln      = log(r1/r2);             % log of their ratio
    if ln == 0, ln = eps; end         % (prevents division by zero)

    % heliocentric angle(s) between target and departure
    acsdp = acos(min(1, max(-1, dotprod)));
    Psi = acsdp; % short way (default)        
    if (longway < 0) , Psi = 2*pi - Psi;  tf = -tf; end    
    if (leftbranch < 0), N = -N; end

    % more initial values
    r1unit = [1, 0, 0];                          % r1 unit (in new coordinates)
    r1th   = [0, 1, 0];                          % tangential direction
    r2vecn = [r2*cos(acsdp), r2*sin(acsdp), 0];  % rotate r2 to new coordinates        
    r2unit = r2vecn/r2;                          % r2 unit (in new coordinates)
    r2th   = [-r2unit(2), r2unit(1), 0];         % tangential direction
    dth    = Psi + 2*pi*N;                       % total turn angle
    k2dth  = k2*dth;                             % easier notation
    k22    = k2*k2;                              % easier notation
    Delta  = 2*(1 - cos(k2dth))/(k2^4) - ln^2;   % Delta term (for gamma_mM)
    k1lim  = min(500, 1/k22);                    % max. permissible [k1] (prevents overflow)
    
    % Delta is negative -- no solution
    if (Delta <= 0), exitflag = -1; return, end

    % limit search space for gamma
    gam1m = atan( k2/2 * ( -ln*cot(k2dth/2) - sqrt(Delta) ));
    gam1M = atan( k2/2 * ( -ln*cot(k2dth/2) + sqrt(Delta) ));
        
    % calculate [k1] factors at the initial points
    k1sgn1 = (ln+tan(gam1m)*sin(k2dth)/k2) / (1-cos(k2dth));
    k11    = sign(k1sgn1) * sqrt(abs(k1sgn1^2 + tan(gam1m)^2/k22));
    k1sgn2 = (ln+tan(gam1M)*sin(k2dth)/k2) / (1-cos(k2dth));
    k12    = sign(k1sgn2) * sqrt(abs(k1sgn2^2 + tan(gam1M)^2/k22));
    gam_interval = (gam1M-gam1m);
    
    % if |k11| > k1lim, adjust gam1m
    if (abs(k11) > k1lim+tol2)
        k1 = inf;  gam1m = gam1m + gam_interval/5;
        addsubtract = sign(k11);   iterations = 0;
        while abs(k1 - addsubtract*k1lim) > tol2
            iterations = iterations + 1;
            k1sgn1     = (ln+tan(gam1m)*sin(k2dth)/k2) / (1-cos(k2dth));
            k1         = sign(k1sgn1)* sqrt(abs(k1sgn1^2 + tan(gam1m)^2/k22));
            dk1dg      = tan(gam1m)*sec(gam1m)^2/k1/k22*(sin(k2dth)/...
                        (1-cos(k2dth))^2*(k2*ln/tan(gam1m)+sin(k2dth))+1);
            gam1m      = gam1m - (k1 - addsubtract*k1lim)/dk1dg;    
            % this function might never fall below the limit; no solution
            if (iterations > 10), exitflag = -1; return, end            
        end
    end
    
    % if |k12| > k1lim, adjust gam1M
    if (abs(k12)> k1lim+tol2)
        k1 = inf;  gam1M = gam1M - gam_interval/5;
        addsubtract = sign(k12);   iterations = 0;
        while abs(k1 - addsubtract*k1lim) > tol2
            iterations = iterations + 1;
            k1sgn1     = (ln+tan(gam1M)*sin(k2dth)/k2) / (1-cos(k2dth));
            k1         = sign(k1sgn1)* sqrt(abs(k1sgn1^2 + tan(gam1M)^2/k22));
            dk1dg      = tan(gam1M)*sec(gam1M)^2/k1/k22*(sin(k2dth)/...
                        (1-cos(k2dth))^2*(k2*ln/tan(gam1M)+sin(k2dth))+1);
            gam1M      = gam1M - (k1 - addsubtract*k1lim)/dk1dg;
            % this function might never fall below the limit; no solution
            if (iterations > 10), exitflag = -1; return, end            
        end
    end
        
    % evaluate the function to test the initial feasibility of the solution.
    % Also evaluate the derivative to test the nature of the graph
    TOFm  = quadrature(gam1m, 'TOF', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);
    TOFpm = quadrature(gam1m, 'derivative', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);
    TOFM  = quadrature(gam1M, 'TOF', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);
    TOFpM = quadrature(gam1M, 'derivative', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);
         
    % at least ONE of these quantities must be positive
    if (TOFm < 0) && (TOFM < 0), exitflag = -1; return, end    
       
    % if both are positive AND both derivatives have equal sign, 
    % the curve lies entirely above zero; no solution possible    
    if (TOFm > 0)&&(TOFM > 0)&&(TOFpm*TOFpM > 0), exitflag = -1; return, end    
    
    % we might have monotonously increasing, monotonously
    % decreasing or "bathtub"-types of solutions. Bathtubs
    % might give rise to two solutions, so first check which 
    % type we're dealing with
    
    % Set initial values for repeated exponential fit routine
    F1    = TOFm;     F2     = TOFM;
    Fp1   = TOFpm;    Fp2    = TOFpM;
    first = gam1m;    second = gam1M;
    
    % monotonously increasing
    if (TOFpm > 0) && (TOFpM > 0)
        % there might not be a solution
        if (TOFm > 0), exitflag = -1; return, end
        
    % monotonously decreasing
    elseif (TOFpm < 0) && (TOFpM < 0)
        % there might not be a solution
        if (TOFM > 0),  exitflag = -1; return, end
        
    % bathtub. Output solution depends on sign of N
    elseif (TOFpm*TOFpM < 0)
                
        % First use the user's setting to check whether there is a
        % solution; the user has selected the left branch (N > 0) or the
        % right branch (N < 0) to be used as a solution. It might be that 
        % for the selected branch, the highest value for the time of flight 
        % falls below zero. In those cases, simply return [exitflag] = -1:
        
        % right branch is below zero; no solution
        if (leftbranch < 0) && (TOFM < 0), exitflag = -1; return, end
        
        % left branch is below zero; no solution
        if (leftbranch > 0) && (TOFm < 0), exitflag = -1; return, end
        
        % If this is NOT the case, we have to find the minimum of the 
        % bathtub in order to determine whether there are solutions or not.
        % This is most easily accomplished by finding the root of the 
        % derivative with Newton-Raphson:

        % fit a simple quadratic function to find an initial estimate
        a   = (TOFpM-TOFpm)/(gam1M-gam1m)/2;
        b   = TOFpm - 2*a*gam1m;
        gmN = -b/2/a; % use the bottom of the parabola as initial estimate
        gmN = min(gam1M, gmN);
        gmN = max(gam1m, gmN);
                
        % Newton-Raphson (simply use forward differences)
        iterations = 0;   h = gam1m*sqrt(eps);   fp = inf;
        while abs(fp) > tol1       
            iterations = iterations + 1;
            fp = quadrature(gmN, 'derivative', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);
            fpp = (quadrature(gmN+h, 'derivative', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]) - fp)/h;
            % naturally, this might fail
            if (fpp == 0), exitflag = -4; return, end
            % Newton-Raphson step
            gmN = gmN - fp/fpp;
            % naturally, this too might fail
            if (iterations > 15) || (gmN < gam1m) || (gmN > gam1M)
                exitflag = -4; return
            end
        end
        
        % evaluate function at (gmN)
        minTOF = quadrature(gmN, 'TOF', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);
         
        % if this value is positive, there are NO solutions
        if (minTOF > 0), exitflag = -1; return, end
        
        % otherwise, select slightly offset initial values for the repeated
        % exponential fit routine (offset is requires because the derivative
        % is of course zero at the minimum. The offset should not be too
        % large though; otherwise, times of flight close to the minimum
        % will not be found):
        
        % right branch
        if (leftbranch < 0)
            first = gmN + (gam1M - gmN)/50;
            F1  = quadrature(first, 'TOF', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);
            Fp1 = quadrature(first, 'derivative', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);
        % left branch
        else
            % first and second points need to be swapped:
            second =  first;  F2 = F1;  Fp2 = Fp1;
            % now evaluate integral at first point           
            first = gmN - (gam1M - gmN)/50;
            F1  = quadrature(first, 'TOF', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);
            Fp1 = quadrature(first, 'derivative', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);            
        end
        
    end
    
    % Use Newton-Raphson in case N = 0
    if (N == 0)
        
        % initial value for Newton-Raphson follows from simple linear fit
        % y = m*x + b
        m  = (TOFM - TOFm)/(gam1M-gam1m);
        b  = TOFm - m*gam1m;        
        gm = -b/m; 
        
        % find gamma with Newton-Raphson
        fx = inf; iterations = 0; 
        while (abs(fx) > tol1)
            % increment iterations
            iterations = iterations + 1; 
            % evaluate function
            fx = quadrature(gm, 'TOF', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);
            fp = quadrature(gm, 'derivative', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);
            % Newton-Raphson
            gm = gm - fx/fp;
            % The straight line fails often; 
            % try exponential-fit: TOF(gamma) = A*exp(B*x) + C;
            if ((gm > gam1M) || (gm < gam1m)) && iterations < 5;
                % make it fail-safe
                if ((Fp1/Fp2) <= 0) || (first == second), exitflag = -4; return, end
                B  = log(Fp1/Fp2)/(first - second);
                A  = (F1 - F2)/(exp(B*first) - exp(B*second));
                C  = min( F1 - A*exp(B*first), F2 - A*exp(B*second));
                % make it fail-safe
                if ((-C/A) <= 0), exitflag = -4; return, end
                % new root
                gm = log(-C/A)/B;
            end
            % the calculation might fail for several reasons
            if (iterations > 25), exitflag = -4; return; end            
        end        
        
    % Use repeated exponential-fit otherwise
    % (it's a lot more stable)
    else
        
        % intialize
        iterations = 0;
        
        % exponential fit method
        while (abs(F1) > tol1) && (abs(F2) > tol1)
            % increment iteration
            iterations = iterations + 1;
                        
            % compute new fit: g(x) = A*exp(B*x) + C
            
            % make it fail-safe
            if ((Fp1/Fp2) <= 0) || (first == second), exitflag = -4; return, end
            B = log(Fp1/Fp2)/(first - second); 
            A = (F1 - F2)/(exp(B*first) - exp(B*second));
            C = min( F1 - A*exp(B*first), F2 - A*exp(B*second));
            if (iterations > 1), first = second; F1 = F2; Fp1 = Fp2; end
            % make it fail-safe
            if ((-C/A) <= 0), exitflag = -4; return, end 
            % compute new root
            second = log(-C/A)/B;
            % evaluate function
            F2  = quadrature(second, 'TOF', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);
            Fp2 = quadrature(second, 'derivative', dth,tf,r1,k2,k2dth,ln,muC,[],[],[],[],[]);
            % This whole method might fail. Return worst exitflag
            if (second > gam1M) || (second < gam1m) || (iterations > 15) 
                exitflag = -4; return;
            end
        end
        
        % copy solution
        gm = second;
        
    end
    
    % get Euler angles for transformation
    crossprd = [r1vec(2)*r2vec(3) - r1vec(3)*r2vec(2),... 
                r1vec(3)*r2vec(1) - r1vec(1)*r2vec(3),...
                r1vec(1)*r2vec(2) - r1vec(2)*r2vec(1)];
    X = r1vec/r1;                                     % direction of X-axis
    Z = crossprd/sqrt(crossprd*crossprd.');           % direction of Z-axis
    beta  = atan2(sqrt(Z(1)^2 + Z(2)^2), Z(3));       % Euler angles
    alpha = atan2(Z(1), -Z(2));
    gamma = atan2(X(3), Z(1)*X(2)-Z(2)*X(1));
    
    % determine rotation matrix for coordinate transformation
    csa = cos(alpha);  csb = cos(beta );   csg = cos(gamma); 
    sna = sin(alpha);  snb = sin(beta );   sng = sin(gamma);
    R = [csa*csg - sng*csb*sna,  csa*csb*sng + sna*csg, snb*sng;
        -sna*csg*csb - csa*sng,  csa*csg*csb - sna*sng, snb*csg;
         snb*sna              , -csa*snb              , csb    ];
    
    % get constants for converged gamma
    [dontcare, k0, k1, phi] = tff(0, gm, r1,k2,k2dth,ln);%#ok
    
    % some frequently used quantities
    k22    = k2^2;       k1k22   = k1*k22;
    k12k24 = k1k22^2;    tangm12 = tan(gm)^2;
    
    % constraint for physically realistic solutions
    % (in principle, this constraint has already been checked before.
    % It's included here, because for some awkward cases this constraint
    % can still be slightly violated after solving the ExpoSin)
    cnstr1 = abs(k1k22) - 1 < tol2;
    
    % constraint that follows from tangential thrust
    % (cos^2 + sin^2 = 1)
    cnstr2  = abs(k12k24 - k22*tangm12 - k12k24*sin(phi)^2) < tol2;       % departure
    tangm22 = (tan(gam1m) + tan(gam1M) - tan(gm))^2;
    cnstr3  = abs(k12k24 - k22*tangm22 - k12k24*sin(k2dth + phi)^2) < tol2;% target
    
    % constraint for k2
    cnstr4 = k22 <= (tangm12+2*k1k22*sin(phi)*ln)/ln^2;
    
    % only if ALL constraints are met, output results
    if cnstr1 && cnstr2 && cnstr3 && cnstr4
        
        % EML DOES NOT SUPPORT ANONYMOUS FUNCTIONS:
        r_start     = k0*exp(k1*sin(phi));
        thdot_start = sqrt(abs(    (muC/r_start^3) /...
                        (tangm12 + k1*k2^2*sin(phi) + 1) ));
        rdot_start  = r_start*thdot_start*k1*k2*cos(phi);
        
        r_end       =  k0*exp(k1*sin(k2*dth + phi));
        thdot_end   =  sqrt(abs(            (muC/r_end^3) /...
                           (tangm22 + k1*k2^2*sin(k2*dth + phi) + 1) ));
        rdot_end    =  r_end*thdot_end*k1*k2*cos(k2*dth + phi);
        
        % components of terminal velocities
        Vth1 = longway * r_start*thdot_start;
        Vr1  = rdot_start;
        Vth2 = longway * r_end*thdot_end;
        Vr2  = rdot_end;
                
        % terminal velocities, transformed back to original space
        V1 = (Vr1*r1unit + Vth1*r1th)*R;
        V2 = (Vr2*r2unit + Vth2*r2th)*R;
        
        % there might be numerical under- or overflow here
        % shouldn't ever happen, but I still encountered it a few times...)
        if any(isnan(V1) | isnan(V2)), exitflag = -3; return, end
        
        % output exposin
        exposin(1) = k0;  exposin( 6) = N;
        exposin(2) = k1;  exposin( 7) = longway*dth;  
        exposin(3) = k2;  exposin( 8) = gm;
        exposin(4) = phi; exposin( 9) = gam1m;                  
        exposin(5) = tf;  exposin(10) = gam1M;
        
        % if we're here, all went fine
        exitflag = 1;
        
        % compute mass
        endmass = quadrature(second, 'endmass', dth,tf,r1,k2,k2dth,ln,muC,k0,k1,phi,M0,Isp);
        
        % compute extremal distance
        extremal_distances = minmax_distances(k0, k1, k2, phi, dth);
        
    else        
        % assign exitflag
        exitflag = -2;        
    end   
    
end % lambert low exposins


% equations for the time of flight
function [tof, k0, k1, phi] = tff(th, gam, r1,k2,k2dth,ln)
          
    % some substitutions
    tangam = tan(gam);
    k22    = k2*k2;
    
    % compute all the constants
    k1sgn = (ln+tangam*sin(k2dth)/k2) / (1-cos(k2dth));
    k1    = sign(k1sgn) * sqrt(abs( k1sgn^2 + tangam^2/k22 ));    
    phi   = acos( min(1, max(-1,tangam/k1/k2)) );
    k0    = r1 * exp(-k1*sin(phi));
        
    % more substitutions
    k2thpphi = k2*th + phi;  
    s        = sin(k2thpphi);           c     = cos(k2thpphi);          
    k12sf    = k1*k22*s;                tangm = k1*k2*c;                
    tanfc    = tangm.^2 + k12sf + 1;    r32   = k0^(1.5)*exp(1.5*k1*s);
        
    % value of the integrand
    tof = r32.*sqrt(abs(tanfc)); % abs is just to be sure
    
end

% equations for the derivative of the time of flight
function tofp = tfp(th, gam, r1,k2,k2dth,ln)
    
    % some substitutions
    tangam = tan(gam);
    k22    = k2^2;
    
    % compute all the constants
    k1sgn = (ln+tangam*sin(k2dth)/k2) / (1-cos(k2dth));
    k1    = sign(k1sgn)* sqrt(abs(k1sgn^2 + tangam^2/k22));
    phi   = real(acos(tangam/k1/k2));
    k0    = r1 * exp(-k1*sin(phi));
    
    % more substitutions
    k2thpphi = k2*th + phi;            cosk2dth = cos(k2dth);
    s        = sin(k2thpphi);          
    c        = cos(k2thpphi);          r32      = k0^(3/2)*exp(3*k1*s/2);
    k12sf    = k1*k22*s;               sec2gam  = 1./cos(gam)^2;
    tangm    = k1*k2*c;                tangam2  = tangam^2;
    tanfc    = tangm.^2 + k12sf + 1;   sink2dth = sin(k2dth);
    
    % prevent one-over zero warnings etc.
    if (k0 == 0), k0 = eps; end
    if any(tanfc == 0)
        % EML DOESN'T SUPPORT LOGICAL INDEXING
        for ii = 1:numel(tanfc)
            if tanfc(ii) == 0, tanfc(ii) = eps; end
        end
    end
        
    % compute the derivatives
    % (yes, they're nasty! But I checked them with central differences,
    % and they're correct)
    dk1dg   = tangam*sec2gam/k1/k22*(sink2dth/(1-cosk2dth)^2*(k2*ln/tangam+sink2dth)+1);
    dphidg  = (dk1dg*tangam-k1*sec2gam)/abs(k1)/sqrt(abs(k22*k1^2-tangam2));
    dk0dg   = -k0*(dk1dg*sin(phi)+dphidg*k1*cos(phi));
    Phi     = dk1dg*s+dphidg*k1*c;
    dPhidth = k2*(dk1dg*c-dphidg*k1*s);
    
    % equation for derivative of transfer time
    tofp = r32.*sqrt(abs(1./tanfc)) .* ( 3*tanfc.*(dk0dg/k0+Phi) + k22*(2*k1*c.*dPhidth/k2+Phi) );
        
end

% · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · 
% Helper functions
% · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · 

% equivalent of Tsjiolkovsky's equation, for low-thrust
% (This gives a mass-estimate after an ExpoSin-leg)
function int = Tsjiolkovsky_low(th, k0, k1, k2, phi)
% Last edited 17/Feb/2009.

    % some substitutions to make life easier
    k22 = k2^2;                      k2thphi = k2*th + phi;
    tangam = k1*k2*cos(k2thphi);     s = sin(k2thphi);
    Gm  = tangam.^2 + k1*k22*s + 1;  r = k0*exp(k1*s);
    cosgam = cos(atan(tangam));

    % the mass integrand    
    int = abs( tangam./cosgam.*(1 - k22*(1 - 2*k1*s)./Gm).*sqrt(1./r./Gm) );
    
end

% Compute minimum and maximum distance to the central body
function extremal_distances = minmax_distances(k0, k1, k2, phi, dth)
        
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
    
    % output distances
    extremal_distances = [minimum_distance, maximum_distance];
    
end

% Perform quadrature with Chebyshev-polynomials; much faster than either
% Gaussian or Newton-Cotes. 
% NOTE: Its accuracy depends on the number of nodes used; this is now 
% determined manually, but in some nearby future this is to be determined 
% automatically. Room for improvement! ^_^
function integral = quadrature(gam, type, dth,tf,r1,k2,k2dth,ln,muC,k0,k1,phi,M0,Isp)

    % Based on the routine given in Numerical Recipes (3rd) section 5.8;
    % calculates the Chebyshev coefficients necessary to approximate some
    % function over the interval [a,b].
    
    % [a,b] = [0,dth]:  (b-a)/2 = (b+a)/2 = dth/2 
    bma = dth/2; bpa = bma;
     
    % select type of integration
    switch type
        case 'TOF'              
            % 50 points suffice
            n = 50; c = zeros(1,n);  k=1:n;
            % evaluate function on Chebychev nodes (note the minus sign; C++ is 0-based)
            y = cos(pi*(k-1/2)/n);   f = tff((y*bma)+bpa, gam, r1,k2,k2dth,ln);
        case 'derivative'
            % 75 points are needed; derivative is much more sensitive
            n = 75;   c = zeros(1,n);  k=1:n;
            % evaluate function on Chebychev nodes (note the minus sign; C++ is 0-based)
            y = cos(pi*(k-1/2)/n);   f = tfp((y*bma)+bpa, gam, r1,k2,k2dth,ln);
        case 'endmass'
            % 50 points suffice also for the mass integral
            n = 50;  c = zeros(1,n);  k=1:n;
            % evaluate function on Chebychev nodes (note the minus sign; C++ is 0-based)
            y = cos(pi*(k-1/2)/n);   f = Tsjiolkovsky_low((y*bma)+bpa, k0,k1,k2,phi);
    end
    
    % compute the coefficients (note the minus sign; C++ is 0-based)
    for j=1:n, c(j)=(f*(cos((pi*(j-1))*((k.'-0.5)/n))))*(2-(j==1))/n; end   

    % An approximation to the integral of the function is
    %   \sum_{k=0}^{n-1}  (c_k · ( k·T_{k+1}(x)/(k^2-1) - x·T_k(x)/(k-1) )) - c(1)·x/2
    % with T_{k}(x) the k-th Chebyshev polynomial of the first kind. Since
    % we need the integral over the interval [-1,1], all the T_k's disappear; 
    %   T_k(1) = 1   and   T_k(-1) = +-1
    % This greatly facilitates our calculation:
    
    F = -[c(1),c(1)]/2;
    for j = 0:(n-1)        
        if (j == 1), F = F + c(2)/2;        
        else            
            F(1) = F(1) + c(j+1)*cos(j*pi)*( 1/(j-1) - j/(j*j-1) ); % x = -1:
            F(2) = F(2) + c(j+1)*( j/(j*j-1) - 1/(j-1) );           % x = +1:
        end
    end
    
    % of course, this should be scaled back from [-1,1] to [a,b]=[0,dth]:
    integral = (F(2)-F(1))*bma; 
    
    % don't forget the constants
    switch type
        case 'TOF'
            integral = integral/sqrt(muC) - tf;
        case 'derivative'
            integral = integral/sqrt(muC)/2;
        case 'endmass' 
            ceff = 9.80665*Isp/1000; % NOTE: km/s             
            integral = M0*exp(-integral*sqrt(muC)/2/ceff);
    end
end
