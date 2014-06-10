function T = ae2T(a, e, muC)
    
    
    
    % get mean motions
    n = ae2n(a, e, muC);
    
    % output periods
    T = zeros(size(e));
    T(e< 1) = 2*pi./n(e<1);
    T(e>=1) = NaN;
    
end
