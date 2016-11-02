function t = aM2t(a, M, t0, muC)
% Author:
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com

% If you find this work useful, please consider a small donation:
% https://www.paypal.me/RodyO/3.5

    % calculate the mean motion
    n = sqrt(muC/abs(a).^3);
    % the times are
    t = M./n - t0;

end