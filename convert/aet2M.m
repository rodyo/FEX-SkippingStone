function M = aet2M(a, e, t, t0, muC)
    
   % get mean motions
   n = ae2n(a, e, muC);
   
   % output is the same for ALL conics
   M = n.*(t - t0)*86400;
   
end