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
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 11/Nov/2009

function E = eM2E(e, M, tol)

    % error check
    error(nargchk(2,3,nargin));%#ok
    if ~all(size(e) == size(M))
        error('eM2E:size_mismatch',...
            'Size of input arguments [e] and [M] must be equal.');
    end

    % set default tolerance
    if (nargin == 2) || isempty(tol), tol = 1e-12; end

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
        
        % apprioximate root from Mikkola's paper
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
