function a = ne2a(n, e, muC)
% NE2A          Convert mean motion and eccentricity to semi-major axis 
%
% a = NE2A(n, e, muC) converts the given mean motion [n] and eccentricity
% [e] to the corresponding semi-major axis [a]. NE2A works correctly for all 
% types of conic sections, and handles scalar/vector/matrix input 
% intuitively. 
%
% NOTE: In case one of the inputs is a parabolic escape trajectory (e = 1), 
% the value taken for [a] is the pericenter distance [rp]. This behavior
% allows all the other conversion routines to work correctly for parabolic
% trajectories.
%
% See also ntheta2t, nM2t, aetheta2t, ae2n.

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
   a   = zeros(size(e)); % initialize output argument
   
   % elliptic cases
   if any(ell(:))
       a(ell) = +1 * (muC./n(ell).^2).^(1/3);
   end
   
   % parabolic cases
   if any(par(:))
       % NOTE: [a] is interpreted as pericenter distance [rp]
       a(par) = (muC./8/n(par).^2).^(1/3);
   end
   
   % hyperbolic cases
   if any(hyp(:))
       a(hyp) = -1 * (muC./n(hyp).^2).^(1/3);
   end    
    
end