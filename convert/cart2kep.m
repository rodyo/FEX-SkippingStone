function varargout = cart2kep(varargin)
% CART2KEP            Convert Cartesian coordinates to Kepler elements
%
% USAGE:
%    [a, e, i, Omega, omega, theta or M] = CART2KEP(statevec, muC, ...
%                                                   'theta' or 'M')
%    [a, e, i, Omega, omega, theta or M] = CART2KEP(x, y, z, dxdt, dydt, ...
%                                              dzdt, muC, 'theta' or 'M')
% or equivalently,
%    elements = CART2KEP(statevec, muC, 'theta' or 'M')
%    elements = CART2KEP(x, y, z, dxdt, dydt, dzdt, muC, 'theta' or 'M')
%
% INPUT ARGUMENTS:
% ================
%   x,    y,    z,   - The Cartesian coordinates to convert. These may
%   dxdt, dydt, dzdt   be arrays of any size, as long as the dimensions
%                      agree.
%           statevec - The following array: [x,y,z, dxdt,dydt,dzdt].
%                muC - The standard gravitational parameter of the central
%                      body.
%     'M' or 'theta' - string argument which decides whether to return the
%                      mean anomaly [M] or true anomaly [theta].
%
% OUTPUT ARGUMENTS:
% =================
%  a, e, i, Omega,   - The kepler elements corresponding to the Cartesian
%  omega, theta or M   coordinates given. Each output has the same size as
%                      the input aruments.
%           elements - The following array: [a,e,i, Omega,omega,theta].
%                      This will be returned if you call CART2KEP() with a
%                      single output argument.
%
%   [a, e, i, O, o, theta] = CART2KEP(x, y, z, xdot, ydot, zdot, muC) or
%   equivalently, CART2KEP( [statevec], muC) will convert the given Cartesian
%   coordinates to Kepler elements. The input value [muC] is a scalar value;
%   the standard gravitational parameter of the central body in question.
%
%   CART2KEP() works correctly for all types of conic sections (circle,
%   ellipse, parabola, hyperbola, rectilinear). The behavior in case of
%   parabolas, zero-inclination orbits and rectilinear paths is somewhat
%   different than regular cases. In case the trajectory turns out to be
%   parabolic, the returned value for [a] is not [inf] (as it is by
%   definition), but the pericenter distance [rp]; this behavior enables
%   all the other conversion routines to handle parabolic cases correctly.
%   For rectilinear paths, the semi-major axis [a] will be zero and the
%   eccentricity [e] will be [inf]. The angles [theta or M], [omega] and
%   [Omega] are ill-defined for rectilinear paths, and they will be
%   returned as [NaN]. The inclination however is still well defined as
%   the angle between the xy-plane and the radius vector (in case V = 0)
%   or velocity vector (in case r = 0). Finally, the angles [omega] and
%   [Omega] are also ill-defined for orbits/trajectories with zero
%   inclination. In such cases, they will simply be taken equal to zero.
%
%   The above is equal to calling CART2KEP(x, y, z, xdot, ydot, zdot, muC,
%   'theta'), that is, the true anomaly [theta] is returned by default. If
%   the mean anomaly is required instead, use CART2KEP(x, y, z, xdot, ydot,
%   zdot, muC, 'M').
%
%   In case the statevector is given as a single array, the Cartesian
%   coordinates are taken columnwise from it. In the more general case of
%   providing each coordinate seperately, the coordinates may be matrices of
%   any dimension, as long as the dimensions are equal for each input.
%
%   Note that if you only require the semi-major axis [a], calling
%
%     a = CART2KEP(...)
%
%   will not work; the expected scalar [a] will be equal to the vector
%
%     [a] = [a,e,i,Omega,omega,theta].
%
%   To return only the semi-major axis, use
%
%     [a,~] = CART2KEP(...)     (MATLAB 2009b and later)
%     [a,dummy] = CART2KEP(...) (MATLAB 2009a and earlier)
%
% instead.
%
% See also kep2cart, cart2para, para2kep, kep2para.

% Author:
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com

% If you find this work useful, please consider a small donation:
% https://www.paypal.me/RodyO/3.5


% NOTE: Although convenient functions, CROSS() and DOT() should be avoided
% in CART2KEP(); When evaluated very often, these two functions produce
% unacceptable overhead and loss of speed.


% Last edited 01/Nov/2009

    % default errortrap
    error(nargchk(2, 8, nargin));%#ok

    % parse input
    thorM = 'M';
    narg = nargin;
    if (narg <= 3)
        statevec = varargin{1};
        % errortraps
        if (size(statevec, 2) < 6)
            error('cart2kep:not_enough_coordinates',...
                  'Statevector should have 6 columns.');
        end
        if (size(statevec, 2) > 6)
            warning('cart2kep:too_many_coordinates',...
                   ['Received statevector with %d columns, while only 6 are required.\n',...
                    'Are you sure this is correct?'], size(statevec, 2));
        end
        % assign coordinates
        szx = [size(statevec,1), 1];
        x = statevec(:, 1); xdot = statevec(:, 4);
        y = statevec(:, 2); ydot = statevec(:, 5);
        z = statevec(:, 3); zdot = statevec(:, 6);
        muC = varargin{2};
    end
    if (narg == 3), thorM = varargin{3}; end
    if (narg >= 7)
        szx = size(varargin{1});
        x = varargin{1}(:);  xdot = varargin{4}(:);
        y = varargin{2}(:);  ydot = varargin{5}(:);
        z = varargin{3}(:);  zdot = varargin{6}(:);
        muC = varargin{7};
    end
    if (narg == 8), thorM = varargin{8}; end

    % errortraps
    if ~isscalar(muC)
        error('cart2kep:muC_not_scalar',...
              'The standard gravitational parameter [muC] should be a scalar value.');
    end
    if ~ischar(thorM)
        error('cart2kep:thorM_not_string',...
              'Selecting the mean anomaly or true anomaly must be done via a string argument.');
    end

    % split position from velocity
    rvec = [x, y, z];
    Vvec = [xdot, ydot, zdot];

    % compute constants required for transformation
    r = sqrt(sum(rvec.^2, 2));   % radius
    zero_r = (r == 0);           % save indices to (r = 0)
    r(zero_r) = realmin;         % prevent division-by-zero
    V = sqrt(sum(Vvec.^2, 2));   % speed
    zero_V = (V == 0);           % save indices to (V = 0)
    H = ...                      % angular momentum vector (=cross(rvec,Vvec,2))
        [rvec(:,2).*Vvec(:,3) - rvec(:,3).*Vvec(:,2),...
         rvec(:,3).*Vvec(:,1) - rvec(:,1).*Vvec(:,3),...
         rvec(:,1).*Vvec(:,2) - rvec(:,2).*Vvec(:,1)];

    % and transform
    N = [-H(:,2), H(:,1),zeros(size(H,1),1)]; % normal vector (to [0,0,1] and H)
    rrr = [r, r, r];                          % for vectorization
    evec = ...                                % eccentricity vector (=cross(Vvec,H,2)/muC-rvec./rrr)
        [Vvec(:,2).*H(:,3) - Vvec(:,3).*H(:,2),...
         Vvec(:,3).*H(:,1) - Vvec(:,1).*H(:,3),...
         Vvec(:,1).*H(:,2) - Vvec(:,2).*H(:,1)]/muC - rvec./rrr;
    e = sqrt(sum(evec.*evec, 2));             % eccentricity
    zero_e = (e == 0);                        % save indices to (e = 0)
    e(zero_e) = realmin;                      % prevent division-by-zero
    Hmag = sqrt(sum(H.*H, 2));                % magnitude of H
    Hmag(Hmag == 0) = NaN;                    % handle rectilinear paths
    i = acos( H(:, 3) ./ Hmag);               % (unambiguous) inclination
    a_divisor = 2./r - V.^2/muC;              % handle parabolic cases
    a = zeros(szx);                           % initialize [a]
    a(a_divisor~=0) = 1./nonzeros(a_divisor); % from vis-viva
    a(a_divisor==0) = (Hmag(a_divisor==0).^2)/2/muC; % parabolic trajectories get [a] = [rp]

    % Omega
    Nxy = sqrt(sum(N.*N, 2));                 % magnitude of normal vector
    Nxy(Nxy == 0) = realmin;                  % prevent division-by-zero
    Omega = atan2(H(:, 1), -H(:, 2));         % Omega is atan2(Hx, -Hy)
    Nxy = [Nxy, Nxy, Nxy];  eee = [e, e, e];  % for vectorization

    % omega, theta
    % cross(N./Nxy, evec,2)
    Ncrossevec = [N(:,2).*evec(:,3) - N(:,3).*evec(:,2),...
                  N(:,3).*evec(:,1) - N(:,1).*evec(:,3),...
                  N(:,1).*evec(:,2) - N(:,2).*evec(:,1)]./Nxy;
    pm1   = sign(sum(Ncrossevec .* H,2)); % sign of the acos for omega
    omega = pm1.*real(acos(sum(evec./eee .* N./Nxy,2))); % omega is acos(e·N) (in units)
    % cross(evec, rvec, 2)
    eveccrossrvec = [evec(:,2).*rvec(:,3) - evec(:,3).*rvec(:,2),...
                     evec(:,3).*rvec(:,1) - evec(:,1).*rvec(:,3),...
                     evec(:,1).*rvec(:,2) - evec(:,2).*rvec(:,1)];
    pm2   = sign( sum( eveccrossrvec .* H, 2 ) );  % sign of the acos for theta
    theta = pm2 .* real(acos( sum(rvec./rrr .* evec./eee, 2) )); % theta is acos(r·e) (in units)

    % corrections for circular orbits
    e(zero_e) = 0;

    % corrections for rectilinear paths
    rectilinear = zero_V | zero_r | all(H == 0, 2) | (e == 1 & a == r/2);
    if any(rectilinear(:))
        a(rectilinear) = 0;    e(rectilinear) = inf;
        i(rectilinear) = atan2( Vvec(rectilinear, 3), ...
            sqrt(Vvec(rectilinear, 1).^2 + Vvec(rectilinear, 2).^2));
        i(rectilinear & i == 0) = atan2( rvec(rectilinear, 3), ...
            sqrt(rvec(rectilinear, 1).^2 + rvec(rectilinear, 2).^2));
        omega(rectilinear) = NaN;
        Omega(rectilinear) = NaN;
        theta(rectilinear)    = NaN;
    end

    % possibly convert theta to M (decided by last input argument)
    if strcmpi(thorM, 'M')
        % parabolic cases handled automatically by etheta2M()
        theta = etheta2M(e, theta);
    end

    % generate apropriate output
    if (nargout == 1)
        varargout{1} = [a, e, i, Omega, omega, theta];
    else
        varargout{1} = reshape(a, szx);   varargout{4} = reshape(Omega, szx);
        varargout{2} = reshape(e, szx);   varargout{5} = reshape(omega, szx);
        varargout{3} = reshape(i, szx);   varargout{6} = reshape(theta, szx);
    end

end
