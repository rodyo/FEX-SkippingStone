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

% Author:
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com

% If you find this work useful, please consider a small donation:
% https://www.paypal.me/RodyO/3.5

    % default errortrap
    error(nargchk(2, 8, nargin));

    % parse & check input
    thorM = 'M';
    narg  = nargin;
    if (narg <= 3)
        elms = varargin{1};
        % error traps
        if (size(elms, 2) < 6)
            error('kep2cart:not_enough_elements',...
                  'Not enough orbital elements; 6 are required, but %d were received.', ...
                  size(elms, 2));
        end
        if (size(elms, 2) > 6)
            warning('kep2cart:too_many_elements',...
                   ['Received too many orbital elements; 6 are required, but %d were received.\n', ...
                    'Are you sure this is correct?'], size(elms, 2));
        end
        % extract input
        a   = elms(:, 1);  O   = elms(:, 4);
        e   = elms(:, 2);  o   = elms(:, 5);
        i   = elms(:, 3);  Mth = elms(:, 6);
        muC = varargin{2};
    end
    if (narg == 3), thorM  = varargin{3}; end
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
    if (narg == 8), thorM  = varargin{8}; end
    if (narg < 7) && (narg > 3)
        error('kep2cart:muC_missing', ['Input argument ''muC'' is missing;\n',...
              'In this mode of operation, KEP2CART() requires at least 7 input arguments.']);
    end

    % select theta or M
    if isempty(thorM), thorM = 'M'; end

    % more errortraps
    if (numel(a) ~= numel(e)) || (numel(a) ~= numel(i)) || (numel(a) ~= numel(O)) || ...
       (numel(a) ~= numel(o)) || (numel(a) ~= numel(Mth))
        error('kep2cart:size_mismatch',...
              'The number of elements in all input matrices must be the same.');
    end
    if ~isscalar(muC)
        error('kep2cart:muC_not_scalar',...
              'The standard gravitational parameter of the central body (muC) must be a scalar.');
    end
    if ~ischar(thorM)
        error('kep2cart:thorM_not_string',...
              'Selection of true or mean anomaly must be done via a string argument.');
    end

    % make it work for [n]-dimensional input
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
    if strcmpi(thorM, 'theta')
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
    H( parabolic) = sqrt( 2*a(parabolic)*muC); % = rp * sqrt(2*muC/rp) = rp * Vesc @ pericenter

    % transform
    xi    = r.*costh;      eta   = r.*sinth;
    one   = [l1; m1; n1] .* [xi; xi; xi];
    two   = [l2; m2; n2] .* [eta; eta; eta];
    rvec  = one + two;
    Vvec  = [-l1.*sinth + l2.*(e+costh)
             -m1.*sinth + m2.*(e+costh)
             -n1.*sinth + n2.*(e+costh)] * muC./[H; H; H];

    % assign values and reshape
    rows = size(rvec, 1) / 3;
    x    = reshape(rvec(1:rows)       , sizea);    xdot = reshape(Vvec(1:rows)       , sizea);
    y    = reshape(rvec(rows+1:2*rows), sizea);    ydot = reshape(Vvec(rows+1:2*rows), sizea);
    z    = reshape(rvec(2*rows+1:end) , sizea);    zdot = reshape(Vvec(2*rows+1:end) , sizea);

    % generate apropriate output
    if (nargout == 1)
        varargout{1} = [x,y,z,xdot,ydot,zdot];
    else
        varargout{1} = x;   varargout{4} = xdot;
        varargout{2} = y;   varargout{5} = ydot;
        varargout{3} = z;   varargout{6} = zdot;
    end

end
