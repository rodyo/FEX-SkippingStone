%#eml   
function [tf, grad_tf, a, grad_a, e] = TOF(theta, rp, r, gradients, muC)
% this function determines the semi-major axis, eccentricity and time
% of flight for given values of [theta], [r] and [rp].

    % default output
    tf = inf; 

    % set initial values (neccessary for EML-compilation)    
    grad_tf = inf(2,1);         grad_a = inf(2,1); 
    dadth = inf;                dadrp = inf;
    dedth = inf;                dedrp = inf;    

    % precompute some quantities
      costh = cos(theta);            sinth = sin(theta);
    costhp1 = 1 + costh;              rpsq = rp*rp;
        rsq = r*r;                     rrp = r*rp;
      denom = r*costhp1 - 2*rp;

    % compute semi-major axis
    a = (rrp*costh - rpsq)/denom;

    % and its gradient
    % (checked with central differences)
    if gradients
        dadth  = sinth*(r*rpsq-rsq*rp)/denom^2;
        dadrp  = (costhp1*(rsq*costh-2*rrp)+2*rpsq)/denom^2;
        grad_a = [dadth; dadrp];
    end

     % compute eccentricity
    e = 1 - rp/a;

    % check if e > 0. If this is not the case, set 
    % everything to [inf] and return
    if (e < 0) || ~isfinite(theta) || ~isfinite(rp)
        grad_tf = [inf; inf]; e = inf;
        a = inf;  grad_a = [inf; inf]; return;
    end

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
    elliptic   = e < 1;  % elliptic case
    parabolic  = e == 1; % parabolic case
    hyperbolic = e > 1;  % hyperbolic case    
    

    % elliptic case
    % (derivatives all checked with central differences)
    if elliptic
        % compute mean motion 
        n = sqrt(muC/a^3);
        % compute the time of flight
        E    = 2*atan2(sqrt1meo1pe*sin(theta/2), cos(theta/2));
        sinE = sin(E);   cosE = cos(E);   tantho2 = tan(theta/2);            
        tf   = (E - e*sin(E))/n;
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
    end
    
    % parabolic case
    % (derivatives all checked with central differences)
    if parabolic
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
    end

     % hyperbolic case
     % (derivatives all checked with central differences)
    if hyperbolic
         % compute mean motion 
        n       = sqrt(muC/(-a)^3);            
        % compute the time of flight
        tantho2 = tan(theta/2);
        F       = 2*atanh( sqrt1meo1pe*tantho2 );
        tf      = (e*sinh(F) - F)/n;    
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
            grad_tf = grad_tf /86400;
        end
        % scale to days
        tf = tf / 86400;

    end % selection of different conics

end % TOF
