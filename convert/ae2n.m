function n = ae2n(a, e, muC)
% AE2N          Convert semi-major axis and eccentricity to mean motion
%
% n = AE2N(a, e, muC) converts the given semi-major axis and eccentricity
% to the corresponding mean motion [n]. AE2N works correctly for all types
% of conic sections, and handles scalar/vector/matrix input intuitively. 
%
% Note that for elliptic and hyperbolic cases, the eccentricity is not
% actually required; it is only there to check if any of the inputs are 
% parabolic. In case one of the inputs is a parabolic escape trajectory 
% (e = 1), it is assumed the corresponding value of [a] is the pericenter
% distance [rp]. This allows computation of the parabolic mean anomaly,
% which is defined in terms of the distance [rp].
%
% See also ntheta2t, nM2t, aetheta2t, ne2a.

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 02/Nov/2009 
    
   % initialize
   ell = (e < 1);   % elliptic orbits
   par = (e == 1);  % parabolic escape trajectories
   hyp = (e > 1);   % hyperbolic escape trajectories
   n   = zeros(size(e)); % initialize output argument
   
   % elliptic cases
   if any(ell(:))
       n(ell) = sqrt(muC./a(ell).^3);
   end
   
   % parabolic cases
   if any(par(:))
       % NOTE: [a] is interpreted as pericenter distance [rp]
       n(par) = sqrt(muC./8/a(par).^3);
   end
   
   % hyperbolic cases
   if any(hyp(:))
       n(hyp) = sqrt(muC./abs(a(hyp)).^3);
   end    
    
end