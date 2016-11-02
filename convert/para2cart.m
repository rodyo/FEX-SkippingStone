function varargout = para2cart(varargin)
% [x, y, z, dxdt, dydt, dzdt] = PARA2CART(a, b, c, u, v, phi) converts the parameterized
% ellipse
%
%    [x y z dxdt dydt dzdt] = [c] + a[u]cos(E) + b[v]sin(E)
%
% for a given value of phi into the corresponding cartesian coordinates.
% This function may also be called by PARA2CART(vec), where vec = [a, b,
% c, u, v, E].
%
% See also

% Author:
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com

% If you find this work useful, please consider a small donation:
% https://www.paypal.me/RodyO/3.5

    % default errortrap
    narg = nargin;
    error(nargchk(2, 7, narg));%#ok

    % parse input parameters
    if (narg == 7)
        params = varargin(1:6);
        muC    = varargin{7};
    end
    if (narg == 2)
        params = varargin{1};
        muC    = varargin{2};
        % possible error
        if ~iscell(params)
            error('para2cart:single_arg_mustbe_cell',...
                  'When calling PARA2CART() with a single argument, that argument must be a cell array.');
        end
        if (size(params,2) ~= 6)
            error('para2cart:size_mismatch',...
                  'Single cell array argument must have 6 columns, received %d.', size(params, 2));
        end
    end

    % continue parsing
    a = params{1}(:);  u = params{4};  if (size(u,2)~=3), u = u.'; end
    b = params{2}(:);  v = params{5};  if (size(v,2)~=3), v = v.'; end
    E = params{6}(:);  c = params{3};  if (size(c,2)~=3), c = c.'; end

    % more error traps
    if (size(u,2)~=3||size(v,2)~=3|| size(c,2)~=3)
        error('para2cart:vectors_badly_sized',...
              'The vectors [u], [v] and [c] all must have one dimension equal to 3.');
    end
    if size(a,1)~=size(b,1)||size(a,1)~=size(E,1)||...
       size(a,1)~=size(u,1)||size(a,1)~=size(v,1)||size(a,1)~=size(c,1)
        error('para2cart:size_mismatch',...
              'Arguments have incompatible size.');
    end

    % calculate eccentricities, mean motions
    e = sqrt(abs(1 - b.*b./a./a));
    n = ae2n(a,e,muC);

    % separate different cases
    ell = (e < 1);  % elliptic orbit
    hyp = (e > 1);  % hyperbolic trajectory

    % replicate matrices to conform to others
    a = a(:, [1,1,1]);  b = b(:, [1,1,1]);   e = e(:, [1,1,1]);
    E = E(:, [1,1,1]);  n = n(:, [1,1,1]);

    % initialize output
    coords = NaN(size(a,1),6);

    % elliptic cases
    if any(ell(:))
        % compute trig functions only once
        cosE = cos(E(ell,:));  sinE = sin(E(ell,:));
        % calculate coordinates
        coords(ell, 1:3) = c(ell,:) + a(ell,:).*u(ell,:).*cosE + b(ell,:).*v(ell,:).*sinE;
        % calculate velocities
        Edot = n(ell,:)./(1 - e(ell,:).*cosE);
        coords(ell, 4:6) = Edot.*(b(ell,:).*v(ell,:).*cosE - a(ell,:).*u(ell,:).*sinE);
    end

    % hyperbolic cases
    if any(hyp(:))
        % compute trig functions only once
        coshF = cosh(E(hyp,:));  sinhF = sinh(E(hyp,:));
        % calculate coordinates
        coords(hyp, 1:3) = c(hyp,:) - a(hyp,:).*u(hyp,:).*coshF + b(hyp,:).*v(hyp,:).*sinhF;
        % calculate velocities
        Fdot = n(hyp,:)./(e(hyp,:).*coshF - 1);
        coords(hyp, 4:6) = Fdot.*( b(hyp).*v(hyp,:).*coshF - a(hyp,:).*u(hyp,:).*sinhF);
    end

    % generate apropriate output
    narg = nargout;
    if (narg == 1)
        varargout{1} = coords;
    elseif (narg == 6)
        varargout{1} = coords(:, 1);   varargout{4} = coords(:, 4);
        varargout{2} = coords(:, 2);   varargout{5} = coords(:, 5);
        varargout{3} = coords(:, 3);   varargout{6} = coords(:, 6);
    end

end % para2cart
