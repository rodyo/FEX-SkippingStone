function theta = aet2theta(a, e, t, t0, muC)

% Author:
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com

% If you find this work useful, please consider a small donation:
% https://www.paypal.me/RodyO/3.5

    % first get M's
    M = aet2M(a, e, t, t0, muC);

    % convert to [theta]
    theta = eM2theta(e, M, 1e-12);

end