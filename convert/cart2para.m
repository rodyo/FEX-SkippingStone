function varargout = cart2para(varargin)
% CART2PARA        Convert cartesian coordinates to parametric 
%                  coefficients
%
% USAGE:
%   [a, b, c, u, v, E0] = CART2PARA(x,y,z,dxdt,dydt,dzdt, muC)
%   [a, b, c, u, v, E0] = CART2PARA(coordinates, muC)
%
%   coefficients = CART2PARA(x,y,z,dxdt,dydt,dzdt, muC)
%   coefficients = CART2PARA(coordinates, muC)
%
% INPUT ARGUMENTS:
% x,y,z,dxdt,dydt,dzdt - Cartesian coordinates of the trajectory to be
%                        converted.
%          coordinates - The following array: [x,y,z,dxdt,dydt,dzdt].
%                 muC  - std. grav. parameter of the central body 
%                        [km3 s-2]
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
% CART2PARA() returns all coefficients [a, b, c, b, u, v, E0] such that 
% the conic section cooresponding to the given Cartesian coordinates is 
% represented by the parametric equation
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
% In case the Cartesian coordinates are given in a single matrix, the
% coordinates [x], [y], [z], [dxdt] and [dydt] will be extracted columnwise. 
% In case each element is given individually, the inputs may be 
% [n]-dimensional, as long as the input matrices contain an equal 
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
% See also para2kep, kep2para, para2cart, kep2cart, cart2kep.

    
% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 14/Nov/2009

    % default errortrap 
    narg = nargin;
    error(nargchk(2, 7, narg));%#ok

    % parse input parameters
    if (narg == 7)
        x = varargin{1}(:);        dxdt = varargin{4}(:);
        y = varargin{2}(:);        dydt = varargin{5}(:); 
        z = varargin{3}(:);        dzdt = varargin{6}(:); 
        muC = varargin{7};  
    elseif (narg == 2)
        cartcoords = varargin{1};
        x = cartcoords(:, 1);   dxdt = cartcoords(:, 4);
        y = cartcoords(:, 2);   dydt = cartcoords(:, 5);
        z = cartcoords(:, 3);   dzdt = cartcoords(:, 6);
        muC = varargin{2};  
    end
    
    % convert to Kepler elements
    elements = cart2kep(x,y,z,dxdt,dydt,dzdt,muC,'M');
    
    % convert to parametric representation
    varargout{1:nargout} = kep2para(elements, muC, 'M');
    
end % cart2para