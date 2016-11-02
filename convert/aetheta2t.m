function t = aetheta2t(a, e, theta, t0, muC)

% Author:
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com

% If you find this work useful, please consider a small donation:
% https://www.paypal.me/RodyO/3.5

    % calculate the mean motions
    n = ae2n(a, e, muC);
    % calculate M
    M = etheta2M(e, theta);
    % the times are given in days
    t = M./n/86400 + t0;

end