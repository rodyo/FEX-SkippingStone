function E = etheta2E(e, th)
% ETHETA2E          Convert eccentricities and true anomalies to 
%                   corresponding eccentric anomalies
%
% ETHETA2E(e, theta) converts the Keplerian elements [e] and [theta] to
% their corresponding eccentric anomaly [E]. The method works correctly for
% all types of conic sections, and handles scalar/vector/matrix input
% intuitively.
% 
% For elliptic orbits, the returned value [E] has the same amount of 
% multiples of 2pi added to it as the original true anomaly [theta]. For
% hyperbolic escape trajectories, the value for [theta] should lie in 
% -acos(-1/e) <= theta <= +acos(-1/e), as can be derived from the physics 
% of such trajectories. If this condition does not hold for the given value
% of [theta], the returned value is simply [NaN]. Finally, for parabolic 
% escape trajectories, there is no analog for the eccentric anomaly [E]. 
% For such cases, the returned value is simply the same true anomaly as is 
% given in the input. 
%
% See also eE2theta, eM2theta.

% Author
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 01/Nov/2009   

    % initialize 
    E = NaN(size(e)); % initialize output argument
    ell = (e < 1);      % elliptic orbits
    par = (e == 1);     % parabolic escape trajectories
    hyp = (e > 1);      % hyperbolic escape trajectories        
    sqrtfac = sqrt(abs(1-e)./(1+e)); % pre-compute the square-root factors

    % elliptic cases
    if any(ell(:))        
        % compute [E] 
        th2pi = mod(abs(th(ell)), 2*pi);
        E(ell) = 2*atan2(sqrtfac(ell).*sin(th2pi/2), cos(th2pi/2));
        % now add the correct amount of 2*pi's to [E] 
        E(ell) = sign(th(ell)).*E(ell) + fix(th(ell)/2/pi)*2*pi;
    end
    
    % parabolic cases
    if any(par(:))
        % there is no analog for [E] in the parabolic case;
        % simply return [theta]
        E(par) = th(par);
    end
    
    % hyperbolic cases
    if any(hyp(:))
        % hyperbolic true anomalies *must* lie in 
        % -acos(-1/e) <= theta <= +acos(-1/e)
        good_th = abs(th(hyp)) <= acos(-1./e(hyp));
        hyp(~good_th) = false;
        % compute [F] (sign uniquely determined by ATANH())
        tano2fc = sqrtfac(hyp).*tan(th(hyp)/2);    
        tano2fc(tano2fc > 1) = 1;   tano2fc(tano2fc < -1) = -1;
        E(hyp) = 2*atanh(tano2fc);
    end        
    
end
