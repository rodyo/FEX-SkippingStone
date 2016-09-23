function varargout = kep2cart(varargin) 
% KEP2CART              Convert Kepler elements to Cartesian coordinates
%
% USAGE:
%    [x, y, z, dxdt, dydt, dzdt] = KEP2CART(elements, muC, 'M' or 'theta') 
%    [x, y, z, dxdt, dydt, dzdt] = KEP2CART(a, e, i, Omega, omega, ...
%                                     M or theta, muC, 'M' or 'theta')
% or equivalently,
%    statevec = KEP2CART(elements, muC, 'M' or 'theta') 
%    statevec = KEP2CART(a, e, i, Omega, omega, ...
%                        M or theta, muC, 'M' or 'theta')
% 
% INPUT ARGUMENTS:
% =================
%  a, e, i, Omega,   - The kepler elements to convert. Each of these should  
%  omega, theta or M   have the same number of elements.
%           elements - The following array: [a,e,i,Omega,omega, theta or M].
%                muC - The standard gravitational parameter of the central
%                      body.     
%     'M' or 'theta' - string argument which decides whether to use the
%                      mean anomaly [M] or true anomaly [theta] for the
%                      conversion.
%
% OUTPUT ARGUMENTS:
% ================
%   x,    y,    z,   - The corresponding Cartesian coordinates. 
%   dxdt, dydt, dzdt   
%           statevec - The following array: [x,y,z, dxdt,dydt,dzdt]. This 
%                      will be returned if you call KEP2CART() with a 
%                      single output argument. 
%
%   [x, y, z, xdot, ydot, zdot] = KEP2CART(a, e, i, O, o, M, muC) or
%   equivalently, KEP2CART( [Keplerian elements], muC) will convert the
%   given orbital elements into Cartesian coordinates. The input value [muC]
%   is a scalar value -- the standard gravitational parameter of the
%   central body in question. This is the same as calling KEPLER2CART(a, e,
%   i, O, o, theta, muC, 'M').
%
%   In case the Keplerian elements are given as a single array, the
%   Keplerian coordinates are taken columnwise from it. In the more general
%   case of providing each element seperately, the elements may be matrices
%   of any dimension, as long as the dimensions are equal for each input.
%
%   Note that the Mean anomaly [M] is used for the default conversion.
%   Calling KEP2CART(a, e, i, o, O, theta, muC, 'theta') (or its
%   matrix equaivalent) will use the true nomaly [theta] directly. 
%
%   KEP2CART() works correctly for all types of conic sections. Parabolic 
%   escape trajectories however form an exception; the value for [a] is 
%   [inf] in theory, which can not provide enough information to complete 
%   the conversion. To circumvent this, it is assumed that when [e] = 1, 
%   the corresponding value for [a] is actually equal to the pericenter 
%   distance, [rp].
%
%  See also cart2kep, kep2para, cart2para.

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com
% Affiliation: LuxSpace sàrl


% If you find this work useful, please consider a donation:
% https://www.paypal.me/RodyO/3.5


    % default errortrap 
    error(nargchk(2,8,nargin, 'struct'));

    % parse & check input
    thorM = 'M';
    narg  = nargin;
    if (narg <= 3)
        elms = varargin{1};        
        assert(size(elms,2) == 6,...
            'cart2kep:size_mismatch',...
            'Element vector should have 6 columns.');
        
        a   = elms(:,1);  O   = elms(:,4);
        e   = elms(:,2);  o   = elms(:,5);
        i   = elms(:,3);  Mth = elms(:,6);                 
        muC = varargin{2};
    end
    if (narg == 3)
        thorM  = varargin{3}; end
    
    if (narg < 2)
        error('kep2cart:muC_missing', ['Input argument ''muC'' is missing;\n',...
              'In this mode of operation, KEP2CART() requires at least 2 input arguments.']);
    end
    if (narg >= 7)
        a = varargin{1};   O = varargin{4}; 
        e = varargin{2};   o = varargin{5};    
        i = varargin{3};   Mth = varargin{6};
        muC = varargin{7};
    end
    if (narg == 8)
        thorM  = varargin{8}; end
    if (narg < 7) && (narg > 3)
        error('kep2cart:muC_missing', ['Input argument ''muC'' is missing;\n',...
              'In this mode of operation, KEP2CART() requires at least 7 input arguments.']);
    end
    
    if isempty(thorM)
        thorM = 'M'; end   
    
    % Basic assertions
    assert(isequal(numel(a),numel(e),numel(i),numel(O),numel(o),numel(Mth)),...
        'kep2cart:size_mismatch',...
        'The number of elements in all input matrices must be equal.');
    assert(isscalar(muC),...
        'cart2kep:muC_not_scalar',...
        'The standard gravitational parameter [muC] should be a scalar value.');    
    assert(ischar(thorM),...
        'cart2kep:thorM_not_string',...
        'Selecting the mean anomaly or true anomaly must be done via a string argument.');

    % make sure it works correctly for [n]-dimensional input
    sizea = size(a);
    a = a(:);    e   = e(:);
    i = i(:);    O   = O(:);
    o = o(:);    Mth = Mth(:);
    
    % indices to parabolic escape trajectories
    parabolic = (e == 1);

    % compute constants required for transformation
    l1 =  cos(O).*cos(o)-sin(O).*sin(o).*cos(i);   l2 = -cos(O).*sin(o)-sin(O).*cos(o).*cos(i);        
    m1 =  sin(O).*cos(o)+cos(O).*sin(o).*cos(i);   m2 = -sin(O).*sin(o)+cos(O).*cos(o).*cos(i);        
    n1 =  sin(o).*sin(i);                          n2 =  cos(o).*sin(i);        
    
    % convert [M] to [theta] (including parabolic cases)    
    if strcmpi(thorM,'theta')
        theta = Mth;
        r = zeros(size(a));
        r(~parabolic) = a(~parabolic).*(1-e(~parabolic).^2) ./...
                       (1 + e(~parabolic).*cos(theta(~parabolic)));
        r( parabolic) = 2*a(parabolic)./(1 + cos(theta(parabolic)));
    else
        % parabolic cases handled automatically by aeM2rtheta()
        [r, theta] = aeM2rtheta(a, e, Mth);
    end
    
    % compute sine & cosine only once
    costh = cos(theta);    
    sinth = sin(theta);
    
    % angular momentum
    H = zeros(size(a));
    H(~parabolic) = sqrt( muC*a(~parabolic).*(1 - e(~parabolic).^2) );
    H( parabolic) = sqrt( 2*a(parabolic)*muC); % = rp*sqrt(2*muC/rp) = rp*Vesc @ pericenter

    % transform    
    xi   = r.*costh;      
    eta  = r.*sinth;
    rvec = [l1; m1; n1].*[xi;xi;xi] + [l2; m2; n2].*[eta;eta;eta];
    Vvec = [-l1.*sinth + l2.*(e+costh) 
            -m1.*sinth + m2.*(e+costh) 
            -n1.*sinth + n2.*(e+costh)] * muC./[H;H;H];

    % assign values and reshape
    rows = size(rvec,1) / 3;
    x    = reshape(rvec(1:rows)       , sizea);    xdot = reshape(Vvec(1:rows)       , sizea);
    y    = reshape(rvec(rows+1:2*rows), sizea);    ydot = reshape(Vvec(rows+1:2*rows), sizea);
    z    = reshape(rvec(2*rows+1:end) , sizea);    zdot = reshape(Vvec(2*rows+1:end) , sizea);        

    % generate appropriate output
    if (nargout == 1)
        varargout{1} = [x,y,z,xdot,ydot,zdot]; 
    else
        varargout{1} = x;   varargout{4} = xdot;
        varargout{2} = y;   varargout{5} = ydot;
        varargout{3} = z;   varargout{6} = zdot;
    end 
    
end


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
% E-mail     : oldenhuis@gmail.com
% Affiliation: LuxSpace sàrl


% If you find this work useful, please consider a donation:
% https://www.paypal.me/RodyO/3.5

        
    % set default tolerance
    if (nargin == 3) || isempty(tol)
        tol = 1e-12; end
    
    % get eccentric and true anomalies
    E  = eM2E(e, M, tol);
    th = eE2theta(e, E);
    
    % initialize
    not_par = (e ~= 1);   % elliptic orbits or hyperbolic escape trajectories
    par = (e == 1);       % parabolic escape trajectories
    r   = zeros(size(e)); % initialize radii 
    
    % radius vectors for non-parabolic trajectories
    if any(not_par(:))
        r(not_par) = ...
            a(not_par).*(1 - e(not_par).^2) ./ ...
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

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com
% Affiliation: LuxSpace sàrl


% If you find this work useful, please consider a donation:
% https://www.paypal.me/RodyO/3.5


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
        th(hyp) = 2*atan2(sqrtfac(hyp).*sinh(E(hyp)/2), cosh(E(hyp)/2)); end
    
end

function E = eM2E(e, M, tol)
% EM2E        Convert eccentricity and mean anomaly to eccentric anomaly
%
% [E] = EM2E(e, M) converts the eccentricity [e] and mean anomaly [M] to 
% the corresponding eccentric anomaly [E]. The mean anomaly [M] should be 
% given in radians. This conversion is essentially a black-box function 
% that solves Kepler's equation for any arbitrary conic section; it works 
% correctly for all types of conic sections, and handles scalar/vector/
% matrix input intuitively. Note that the eccentric anomaly has no
% definition in the parabolic case; for parabolic cases, [E] = [theta] 
% (the true anomaly) is returned. 
%
% EM2E(e, M, tol) (with a third argument) uses a maximum error-tolerance 
% given by [tol]. The default is 1e-12.
%
% Algorithm: To calculate [theta] from [M], EM2E uses a carefully selected 
% first approximate root of Kepler's equation, followed by a Newton-Raphson 
% iteration scheme if the first estimate is not within the limits set by 
% [tol]. See [Seppo Mikkola, "A cubic approximation for Kepler's equation", 
% 1987] for more details.
%
% See also eM2theta, eE2theta, eE2M.

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com
% Affiliation: LuxSpace sàrl


% If you find this work useful, please consider a donation:
% https://www.paypal.me/RodyO/3.5


    % error check
    error(nargchk(2,3,nargin,'struct'));
    assert(all(size(e) == size(M),...
        'eM2E:size_mismatch',...
        'Size of input arguments [e] and [M] must be equal.');
    
    % set default tolerance
    if (nargin == 2) || isempty(tol)
        tol = 1e-12; end

    % initialize
    sze = size(e);  % save original size of [e]
    M   = M(:);     % convert to colum vector       
    e   = e(:);     % convert to colum vector              
    E   = zeros(size(e)); % initialize output argument   
    ell = (e < 1);  % elliptic orbits     
    par = (e == 1); % parabolic escape trajectories
    hyp = (e > 1);  % hyperbolic escape trajectories
   
    % elliptic cases
    if any(ell)

        % extract appropriate values
        Me = M(ell);  ee = e(ell);

        % put M in interval (-pi < M < +pi)
        if any(Me > pi) || any(Me < -pi)
            Me = mod(Me, 2*pi);
            Me(Me > pi) = Me(Me > pi) - 2*pi;            
        end
        
        % Apprioximate root from Mikkola's paper
        gamma = (4*ee + 1/2);
        alpha = (1 - ee)./gamma;
        beta  = Me./(2*gamma);
        z     = (beta + sign(beta).*sqrt(beta.^2 + alpha.^3));
        z     = sign(z).*abs(z).^(1/3);
        z(z==0) = realmin;  % prevent 1-over-0 
        s  = z - alpha./z;
        s  = s - 0.078*s.^5./(1 + ee);
        Ee = Me + ee.*(3*s - 4*s.^3);
        Ep = inf;

        % perform Newton-Raphson iterations
        while any(abs(Ee-Ep) >= tol)
            Ep = Ee;  
            Ee = Ee + (Me + ee.*sin(Ee) - Ee) ./ (1 - ee.*cos(Ee));
        end

        % make sure special cases are exact
        Ee(Me ==   0) = 0;
        Ee(Me == +pi) = +pi;
        Ee(Me == -pi) = -pi;

        % add the same amount of multiples of 2pi to [E] as 
        % there were in [M]        
        Ee = Ee + sign(M(ell)).*fix((abs(M(ell))+pi)/2/pi)*2*pi;

        % insert into final array
        E(ell) = Ee;

    end

    % hyperbolic cases
    if any(hyp)

        % extract relevant values
        eh = e(hyp); Mh = M(hyp);

        % apprioximate root from Mikkola's paper
        gamma = (4*eh + 1/2);
        alpha = (eh - 1)./gamma;
        beta  = Mh./(2*gamma);
        z  = (beta + sign(beta).*sqrt(beta.^2 + alpha.^3));
        z  = sign(z).*abs(z).^(1/3);
        s  = z - alpha./z;
        s2 = s.^2;
        s  = s + 0.071*s.^5 ./ ((1 + 0.45.*s2).*(1 + 4*s2).*eh);
        Eh = 3*log(s + sqrt(1 + s.^2));

        % perform at least one Newton-Raphson iteration, additional ones if
        % required by [tol]
        Ep = inf;
        while any( abs(Eh - Ep) >= tol)
            Ep = Eh;
            Eh = Eh + (Eh - eh.*sinh(Eh) + Mh) ./ (eh.*cosh(Eh) - 1);
        end

        % insert into final array
        E(hyp) = Eh;            

    end   
    
    % solve parabolic cases
    if any(par(:))
        % [E] is not defined; return [theta]
        E(par) = eM2theta(e(par), M(par));
    end

    % transform back to original size
    E = reshape(E, sze);

end

function th = eM2theta(e, M, tol)  
% EM2THETA      Convert eccentricity and mean anomaly to true anomaly
%
% [th] = EM2THETA(e, M) converts the the eccentricity [e] and mean 
% anomaly [M] to the corresponding true anomaly [th]. The mean anomaly 
% [M] should be given in radians. This conversion is essentially an
% black-box function that solves Kepler's equation for any arbitrary 
% conic section; it works correctly for all types of conic sections, 
% and handles scalar/vector/matrix input.
%
% EM2THETA(e, M, tol) (with a third argument) uses a maximum error 
% tolerance for the iterative procedure given by [tol]. The default 
% is 1e-12.
%
% Algorithm: To calculate [theta] from [M], EM2THETA uses a carefully 
% selected first approximate root of Kepler's equation, followed by a 
% Newton-Raphson iteration scheme if the first estimate is not within 
% the limits set by [tol]. See [Seppo Mikkola, "A cubic approximation 
% for Kepler's equation", 1987] for more details.
%
% See also eM2E, eE2M, eE2theta.


% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com
% Affiliation: LuxSpace sàrl

% If you find this work useful, please consider a donation:
% https://www.paypal.me/RodyO/3.5

    % set default tolerance        
    if (nargin == 2) || isempty(tol), tol = 1e-12; end

    % initialize       
    par = (e == 1);     % parabolic cases        
    th  = NaN(size(e)); % initialize output argument

    % get eccentric anomalies and calculate theta for 
    % non-parabolic cases
    if ~all(par(:)) 
        Enp = eM2E(e(~par), M(~par), tol);
        th(~par) = eE2theta(e(~par), Enp);        
    end

    % handle parabolic cases
    % (analytic solution to Barker's equation)
    if any(par(:))            
        % extract proper [M]'s
        Mp = M(par);     
        % prevent 1-over-zero warnings
        Mp(Mp == 0) = realmin;
        % compute [th]
        y = abs(atan(1/3./Mp));
        x = atan(tan(y/2).^(1/3));
        th(par) = 2*atan(2./tan(2*x));            
        % make sure special cases are exact
        th(par & M==0) = 0;
        % make sure the quadrant is also correct
        th(par) = th(par).*sign(Mp);            
    end
        
end
