function theta = aet2theta(a, e, t, t0, muC)
    
    % first get M's
    M = aet2M(a, e, t, t0, muC);
    
    % convert to [theta]
    theta = eM2theta(e, M, 1e-12);
    
end