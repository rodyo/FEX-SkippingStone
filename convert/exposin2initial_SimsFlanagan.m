% velocity, thrust, position and Delta-V of an ExpoSin, for an
% equally-spaced amount of points [steps], point-to-point.
%
% ( Serves as an initial estimate for SIMSFLANAGAN() ) 
function [burnprog, thetas, ts] = ...
        exposin2initial_SimsFlanagan(exposin, r1vec, r2vec, dt, steps)
% EXPOSIN2INITIAL_DELTAV

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 24/Mar/2009.

    % get global
    global muS

    % times vector
    t = linspace(0, dt, steps);

    % transformation matrix
    X = unit(r1vec, 2);
    Z = unit(cross(r1vec, r2vec, 2), 2);
    Y = cross(Z, X);
    beta  = atan2(sqrt(Z(1)^2 + Z(2)^2), Z(3));
    alpha = atan2(Z(1), -Z(2));
    gamma = atan2(X(3), Y(3));
    R = angle2dcm(alpha, beta, gamma, 'ZXZ');

    % extract exposin coefficients
    k0  = exposin(1);      phi = exposin(4);
    k1  = exposin(2);
    k2  = exposin(3);      dth = exposin(7);

    % substitutions
    tgm   = @(th) k1*k2*cos(k2*th + phi);
    s     = @(th) sin(k2*th + phi);
    d     = @(th) (tgm(th).^2 + k1*k2^2*s(th) + 1);

    % define functions for position and velocity
    r     = @(th) k0.*exp(k1.*s(th));
    thjdot = @(t, th) sqrt( ( muS./r(th).^3) ./ d(th) );
    rdot  = @(th) r(th) .* tgm(th) .* thjdot(0, th);

    % calculate points
    [ts, thj] = ode113(thjdot, t, 0);
    rj        = r(abs(thj));

    % find velocities at all points (j)
    runit  = [cos(thj),  sin(thj), zeros(size(rj))];
    thunit = [-sin(thj), cos(thj), zeros(size(rj))];
    Vr     = rdot(thj);
    Vth    = rj.*thjdot(0, thj);
    Vjvecs = [Vr, Vr, Vr].*runit + [Vth, Vth, Vth].*thunit;
    Vjs    = sqrt(sum(Vjvecs.*Vjvecs, 2));

    % determine DeltaVs
    Vjunits  = unit(Vjvecs, 2);
    energies = Vjs.^2/2 - muS./rj;
    energies = circshift(energies, -1);
    Vnews = sqrt(2*(energies + muS./rj));
    deltaV = Vnews - Vjs;
    deltaV(end) = 0;
    deltaVs = [deltaV, deltaV, deltaV] .* Vjunits;

    % final output
    burnprog = deltaVs*R;
    thetas   = thj;

end
