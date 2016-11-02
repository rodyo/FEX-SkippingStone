function M = aet2M(a, e, t, t0, muC)

% Author:
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com

% If you find this work useful, please consider a small donation:
% https://www.paypal.me/RodyO/3.5

   % get mean motions
   n = ae2n(a, e, muC);

   % output is the same for ALL conics
   M = n.*(t - t0)*86400;

end