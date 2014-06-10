function M = etheta2M(e, theta)
% ETHETA2M          Convert eccentricities and true anomalies to 
%                   corresponding mean anomalies
%
% ETHETA2M(e, theta) converts the Keplerian elements [e] and [theta] to
% their corresponding mean anomalies [M]. The method works correctly for
% all types of conic sections, and handles scalar/vector/matrix input 
% intuitively.
% 
% For elliptic orbits, the returned value [M] has the same amount of 
% multiples of 2pi added to it as the original true anomaly [theta]. For
% hyperbolic escape trajectories, the value for [theta] should lie in 
% -acos(-1/e) <= theta <= +acos(-1/e), as can be derived from the physics 
% of such trajectories. If this condition does not hold for the given value
% of [theta], the returned value is simply [NaN]. Finally, for parabolic 
% escape trajectories, the mean anomaly returned follows from Barker's
% equation;  tan(th/2) + 1/3*tan(th/2)^3 = 2M.

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 01/Nov/2009 

    % seperate ellipses from hyperbolas
    M   = zeros(size(e)); % initialize output argument
    ell = (e < 1);        % elliptic orbits
    par = (e == 1);       % parabolic escape trajectories
    hyp = (e > 1);        % hyperbolic escape trajectories
    
    % first compute eccentric anomalies
    if ~all(par) % ([E] not needed if all given cases are parabolic)
        E = etheta2E(e, theta);
    end
    
    % elliptic cases
    if any(ell(:))
        M(ell) = E(ell) - e(ell).*sin(E(ell));
    end
    
    % parabolic cases (from Barker's equation)
    if any(par(:))
        tantho2 = tan(theta(par)/2);
        M(par)  = (tantho2 + (tantho2.^3)/3)/2;
    end
    
    % hyperbolic cases
    if any(hyp(:))
        M(hyp) = e(hyp).*sinh(E(hyp)) - E(hyp);
    end

end
