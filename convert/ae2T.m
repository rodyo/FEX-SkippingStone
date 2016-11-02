function T = ae2T(a, e, muC)

% Author:
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com

% If you find this work useful, please consider a small donation:
% https://www.paypal.me/RodyO/3.5


    % get mean motions
    n = ae2n(a, e, muC);

    % output periods
    T = zeros(size(e));
    T(e< 1) = 2*pi./n(e<1);
    T(e>=1) = NaN;

end
