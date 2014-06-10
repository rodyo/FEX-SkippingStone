function th = eE2theta(e, E)
% EE2THETA          Convert eccentric anomalies and eccentricities to
%                   the corresponding true anomalies %                    
%
% EE2THETA(e, E) converts the eccentricity [e] and eccentric anomaly [E] to
% the corresponding true anomaly [theta]. The method works correctly for
% all types of conic sections, and handles scalar/vector/matrix input
% intuitively.
% 
% For elliptic orbits, the returned value [th] has the same amount of 
% multiples of 2pi added to it as the original eccentric anomaly [E]. 
% Finally, for parabolic escape trajectories there is no analog for 
% the eccentric anomaly [E]. If [e] = 1 and an [E] is given, the 
% returned value will simply be equal to [E].
%
% See also etheta2E, etheta2M, eM2theta.

% Author
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 01/Nov/2009   

    % initialize
    th = zeros(size(e)); % pre-allocate theta
    ell = (e < 1);       % elliptic orbits
    par = (e == 1);      % parabolic escape trajectories
    hyp = (e > 1);       % hyperbolic escape trajectories
    sqrtfac = sqrt((1+e)./abs(1-e)); %pre-compute squareroot factors

    % elliptic cases
    if any(ell(:))
        % compute thetas
        E2pi = mod(abs(E(ell)), 2*pi);
        th(ell) = 2*atan2(sqrtfac(ell).*sin(E2pi/2), cos(E2pi/2));
        % add correct amount of multiples of 2pi to [theta]
        th(ell) = sign(E(ell)).*th(ell) + fix(E(ell)/2/pi)*2*pi;
    end
    
    % parabolic cases
    if any(par(:))
        % [E] undefined for parabolic cases; return same values
        th(par) = E(par);
    end
    
    % hyperbolic cases
    if any(hyp(:))
        th(hyp) = 2*atan2(sqrtfac(hyp).*sinh(E(hyp)/2), cosh(E(hyp)/2));
    end
    
end
