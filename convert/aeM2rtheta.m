function [r, th] = aeM2rtheta(a, e, M, tol)  
% AEM2RTHETA        Convert semimajor axis, eccentricity and mean anomaly
%                   to corresponding polar coordinates
%
% [r, th] = AEM2RTHETA(a, e, M) converts the semi-major axis [a], the
% eccentricity [e] and mean anomaly [M] to the polar coordinates [r] and
% [th]. The mean anomaly [M] should be given in radians. This conversion 
% is essentially a black-box function that solves Kepler's equation for 
% any arbitrary conic section; it works correctly for all types of conic 
% sections, and handles scalar/vector/matrix input.
%
% AEM2RTHETA(a, e, M, tol) uses a maximum error-tolerance given by [tol].
% The default is 1e-12.
%
% NOTE: for parabolic trajectories, the semi-major axis [a] is infinite and
% does not provide enough information to complete the conversion. For that 
% reason it is assumed that for all [e] == 1, the corresponding value for
% [a] is NOT [inf], but equal to the pericenter distance [rp]; doing so 
% DOES provide enough information to compute [r] and [th].
%
% Algorithm: To calculate [theta] from [M], AEM2RTHETA() uses a carefully 
% selected first approximate root of Kepler's equation, followed by a 
% Newton-Raphson iteration scheme if the first estimate is not within the 
% limits set by [tol]. See [Seppo Mikkola, "A cubic approximation for 
% Kepler's equation", 1987] for more details.
%
% See also eM2theta, eE2theta, etheta2M.

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 02/Nov/2009
        
    % set default tolerance
    if (nargin == 3) || isempty(tol), tol = 1e-12; end
    
    % get eccentric and true anomalies
    E  = eM2E(e, M, tol);
    th = eE2theta(e, E);
    
    % initialize
    not_par = (e ~= 1);   % elliptic orbits or hyperbolic escape trajectories
    par = (e == 1);       % parabolic escape trajectories
    r   = zeros(size(e)); % initialize radii 
    
    % radius vectors for non-parabolic trajectories
    if any(not_par(:))
        r(not_par) = a(not_par).*(1 - e(not_par).^2) ./ ...
                    (1 + e(not_par).*cos(th(not_par)));
    end
    
    % radius vectors for parabolic trajectories
    if any(par(:))
        % assuming [a](par) = [rp], we have [p] = 2*[rp].
        % Then, using cos^2(gamma) = (1 + cos(th))/2, and
        % [p] = 2*r*cos^2(gamma), the following is easily derived:
        r(par) = 2*a(par)./(1 + cos(th(par)));
    end
    
end
