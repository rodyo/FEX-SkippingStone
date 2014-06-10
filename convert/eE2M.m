function M = eE2M(e, E)
% EE2M        Convert eccentricity and eccentric anomaly to mean anomaly
%
% [M] = EE2M(e, E) converts the eccentricity [e] and eccentricity anomaly 
% [E] to the corresponding mean anomaly [M]. The eccentric anomaly [E] 
% should be given in radians. EE2M() works correctly for all types of conic 
% sections, and handles scalar/vector/matrix input intuitively. Note that 
% the eccentric anomaly has no definition in the parabolic case; for 
% parabolic cases, it is assumed that the given [E] is the true anomaly 
% [theta], and the returned value for [M] follows from Barker's equation.
%
% See also eM2theta, eE2theta, eM2E.

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 02/Nov/2009       

    % initialize
    ell = (e < 1);   % elliptic orbits
    par = (e == 1);  % parabolic escape trajectories
    hyp = (e > 1);   % hyperbolic escape trajectories
    M   = zeros(size(e)); % initialize output argument

    % calculate mean anomalies
    if any(ell(:))
        M(ell) = E(ell) - e(ell).*sin(E(ell));
    end
    if any(par(:))
        tantho2 = tan(E(par)/2);
        M(par)  = (tantho2 + (tantho2.^3)/3)/2;
    end
    if any(hyp(:))
        M(hyp) = e(hyp).*sinh(E(hyp)) - E(hyp);
    end
    
end    
