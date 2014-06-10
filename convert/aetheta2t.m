function t = aetheta2t(a, e, theta, t0, muC)
    
    
    % calculate the mean motions
    n = ae2n(a, e, muC);
    % calculate M
    M = etheta2M(e, theta);
    % the times are given in days
    t = M./n/86400 + t0;
    
end