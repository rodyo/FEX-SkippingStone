function statevec = two_body_ephemerides(body, time, t0, model, constants)
%
    
    % extract orbital elements
    a = model.as    (body);  o = model.omegas(body);
    e = model.es    (body);  M = model.M0s   (body);
    i = model.is    (body);  n = model.ns    (body);
    O = model.Omegas(body);
        
    % put M at the desired time
    M = M + n*((time-t0)*86400); 
        % (convert to seconds, as [time] is given in days)
    
    % return the coordinates at this time
    [x, y, z, xdot, ydot, zdot] = kep2cart(a, e, i, O, o, M, constants.mu_central, 'M');
    statevec = [x, y, z, xdot, ydot, zdot];
    
end