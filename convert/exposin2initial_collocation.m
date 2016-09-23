% mass, velocity, thrust and position of an ExpoSin, for an
% equally-spaced amount of points [steps].
%
% ( Serves as an initial estimate for COLLOCATION() ) 
function initial = exposin2initial_collocation(exposin, r1vec, r2vec, m0, mf, steps)
%EXPOSIN2INITIAL

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 24/Mar/2009.



    % extract exposin coefficients
    k0  = exposin(1);      phi = exposin(4);
    k1  = exposin(2);
    k2  = exposin(3);      dth = exposin(7);

    % substitutions
    tgm   = @(th) k1*k2*cos(k2*th + phi);
    s     = @(th) sin(k2*th + phi);
    d     = @(th) (tgm(th).^2 + k1*k2^2*s(th) + 1);

    % define functions for position and velocity
    r      = @(th) k0.*exp(k1.*s(th));
    thjdot = @(t, th) sqrt( ( muS./r(th).^3) ./ d(th) );
    rdot   = @(th) r(th) .* tgm(th) .* thjdot(0, th);

    % define functions for thrust
    aN = @(th) tgm(th) / 2 ./ cos(atan(tgm(th)));
    an = @(th) aN(th) .* (1./d(th) - (k2^2*(1 - 2*k1*s(th)))./d(th).^2);
    at = @(th) an(th) .* (muS ./ r(th).^2);

    % calculate points
    r   = @(th) k0.*exp(k1.*sin(k2*th + phi));
    thj = linspace(0, dth, steps).';
    rj  = r(abs(thj));
    [x, y] = pol2cart(thj, rj);

    % coordinate transformation
    X = unit(r1vec, 2);
    Z = unit(cross(r1vec, r2vec, 2), 2);
    Y = cross(Z, X);
    beta  = atan2(sqrt(Z(1)^2 + Z(2)^2), Z(3));
    alpha = atan2(Z(1), -Z(2));
    gamma = atan2(X(3), Y(3));
    R = angle2dcm(alpha, beta, gamma, 'ZXZ');
    [xj, yj, zj] = blkassign(([x, y, zeros(size(x))]*R));

    % find velocities at all points (j)
    runit   = [ cos(thj), sin(thj), zeros(size(rj))];
    thjunit = [-sin(thj),  cos(thj), zeros(size(rj))];
    Vr      = rdot(abs(thj));
    Vthj    = rj.*thjdot(0, abs(thj));
    Vjvecs  = [Vr, Vr, Vr].*runit + [Vthj, Vthj, Vthj].*thjunit;
    Vj      = Vjvecs*R;

    % find masses at all points (j)
    massj = (linspace(m0, mf, steps)).';
    massj = 1500*ones(size(massj));

    % find thrusts at all points (j)
    acc    = at(thj);
    tgunit = unit(Vj, 2);
    accj   = tgunit .* [acc, acc, acc];

    % final output
    initial = [flipud(xj), flipud(yj), flipud(zj), flipud(-Vj), massj, accj].';
    initial = initial(:);

end
