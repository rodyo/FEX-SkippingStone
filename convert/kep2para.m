function varargout = kep2para(varargin)
% KEP2PARA        Convert kepler elements to parametric representation
%
% USAGE:
%   [a, b, c, u, v, E0] = KEP2PARA(...
%               a,e,i,Omega,omega,theta or M, muC, 'theta' or 'M')
%
%   [a, b, c, u, v, E0] = KEP2PARA(elements, muC, 'theta' or 'M')
%
%   coefficients = KEP2PARAKEP2PARA(...
%               a,e,i,Omega,omega,theta or M, muC, 'theta' or 'M')
%
%   coefficients = KEP2PARA(elements, muC, 'theta' or 'M')
%
%
% INPUT ARGUMENTS:
%   a,e,i,Omega,omega, - Kepler elements of the conic section to be
%          theta or M    converted. Note that either the true anomaly
%                        [theta] or mean anomaly [M] can be provided. All
%                        elements may be given as matrices of size (Nx1).
%             elements - The following array: [a,e,i,Omega,omega,theta/M].
%                 muC  - std. grav. parameter of the central body 
%                        [km3 s-2]
%       'theta' or 'M' - string arument which indicates whether the true
%                        anomaly [theta] or mean anomaly [M] was provided.
%                        Defaults to 'M' (mean anomaly).
%
% OUTPUT ARGUMENTS:
%         a,b,c,u,v,E0 - Coefficients for the parametric representation.
%                        They are:
%                            a: semi-major axis           (Nx1 scalars)
%                            b: semi-minor axis           (Nx1 scalars)
%                            c: location of the center    (Nx3 vectors)
%                            u: unit vector to pericentr  (Nx3 vectors)
%                            v: unit vector to [b]        (Nx3 vectors)
%                           E0: initial eccentric anomaly (Nx1 scalars)
%         coefficients - The following cell-array: {a, b, c, u, v, E0}.
%
% KEP2PARA() returns all coefficients [a, b, c, b, u, v, E0] such that 
% the corresponding conic section is represented by the parametric 
% equation
%
%   ellipses:
%    [x y z dxdt dydt dzdt] = ...
%               [ [c] + a[u]cos(E)  + b[v]sin(E)  (position)
%      n/(1-e*cos(E))*(-a[u]sin(E)  + b[v]cos(E)] (velocity)
%
%   hyperbolae:
%    [x y z dxdt dydt dzdt] = ...
%               [ [c] - a[u]cosh(E) + b[v]sinh(E) (position)
%     n/(e*cosh(E)-1)*(-a[u]sinh(E)  + b[v]cosh(E)] (velocity)
%
% Where [n] is the mean motion and [E] (the eccentric anomaly) the 
% parameter. This representation is sometimes useful when an extremely 
% large amount of statevectors needs to be known for a particular conic
% section; the parametric representation eliminates the need to convert
% the eccentric anomaly to true anomaly before the statevector can be
% evaluated. Moreover, when the conic section has low eccentricity 
% (e <~ 0.3), this parametric representation provides a way to quickly 
% compute an estimate of the statevector at a particular instance, by 
% using the approximation 
%
%    E(t) ~ 2*pi*t/T + E0
%
% in the parametric equation above. Here, [T] is the orbital period, 
% and [t] time. The accuracy of this approximation is exact for [e] = 0, 
% and quickly grows worse with increasing eccentricity. The aforementioned 
% eccentricity (e = 0.3) constitutes to a maximum error of 10% in both 
% position and velocity. 
%
% In case the Keplerian elements are given in a single matrix, the
% elements [a], [e], [i], [omega] and [Omega] will be extracted
% columnwise. In case each element is given individually, the inputs may
% be [n]-dimensional, as long as the input matrices contain an equal 
% amount of elements [N]. In either case, the corresponding output always 
% have size 
%
%   a  (N x 1)      c  (N x 3)  
%   b  (N x 1)      u  (N x 3)
%   E0 (N x 1)      v  (N x 3)
%
% In case a single output is requested (all coefficients in a single output
% argument), the output argument will be a cell-array, containing
%
%   {a, b, c, u, v, E0} = {[Nx1] [Nx1] [Nx3] [Nx3] [Nx3] [Nx1]}
%   
% See also para2kep, cart2para, para2cart, kep2cart, cart2kep.


% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 14/Nov/2009

    % default errortrap 
    narg = nargin;
    error(nargchk(2, 8, narg));%#ok

    % parse input parameters
    thorM = 'M'; % defaults to mean anomaly
    if (narg >= 7)
        a = varargin{1}(:);        O    = varargin{4}(:);
        e = varargin{2}(:);        o    = varargin{5}(:); 
        i = varargin{3}(:);        Mth0 = varargin{6}(:); 
        muC = varargin{7};  
        if (narg == 8), thorM = varargin{8}; end
    elseif (narg <= 3)
        elements = varargin{1};
        a = elements(:, 1);   O    = elements(:, 4);
        e = elements(:, 2);   o    = elements(:, 5);
        i = elements(:, 3);   Mth0 = elements(:, 6);
        muC = varargin{2};
        if (narg == 3), thorM = varargin{3}; end
    end

    % determine semi-minor axes
    b = a .* sqrt(abs(1 - e.^2)); % both ellipses and hyperbolas

    % determine rp, ra
    [rpx, rpy, rpz] = kep2cart(a, e, i, O, o, zeros(size(a)), muC, 'theta'); % thorM not needed
    rp = a.*(1 - e); 
    rpunit = bsxfun(@rdivide, [rpx, rpy, rpz], rp);

    % coordinates of center depends on shape of the conic;
    % if it's elliptic, the product a*e should be subtracted from
    % pericenter. If it's hyperbolic, it should be added:
    c = bsxfun(@times, sign(e-1).*a.*e, rpunit);

    % calculate [c], [u], [v], [w]
    u = rpunit;
    [vx, vy, vz] = kep2cart(a, e, i, o, O, pi/2-e, muC, 'theta');
    v = [vx, vy, vz] - c;
    v = bsxfun(@rdivide, v, sqrt(sum(v.^2, 2)));
    
    % calculate [E0]
    if strcmpi(thorM, 'M')
        E0 = eM2E(e, Mth0);
    elseif strcmpi(thorM, 'theta')
        E0 = etheta2E(e, Mth0);
    end

    % generate apropriate output
    narg = nargout;
    if (narg == 1)
        varargout{1} = {a, b, c, u, v, E0};
    elseif (narg > 1)
        varargout{1} = a;   varargout{4} = u;
        varargout{2} = b;   varargout{5} = v;
        varargout{3} = c;   varargout{6} = E0;
    end

end % kep2para
