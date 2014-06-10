function t = aM2t(a, M, t0, muC)
    
    
    % calculate the mean motion
    n = sqrt(muC/abs(a).^3);
    % the times are
    t = M./n - t0;
    
end