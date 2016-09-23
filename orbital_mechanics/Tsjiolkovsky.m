function out = Tsjiolkovsky(DeltaV, Isp, M0, Me)
%TSJIOLKOVSKY       Determine parameters from Tsjiolkovsky's equation.
%
%   value = TSJIOLKOVSKY(DeltaV, Isp, M0, Me) should be issued with one
%   of the input parameters set to -1. The function then calculates its
%   value. TSJIOLKOVSKY_HIGH is basically an implementation of all possible
%   ways to solve Tsjiolkovsky's equation, once the three other
%   parameters are known.
%
%   TSJIOLKOVSKY_HIGH is vectorized in the sense that all inputs may be arrays,
%   with one input equal to a scalar value of -1. In that case, all
%   corresponding values for the requested output will be calculated.
%
%   Units are assumed to be (kg) for the masses, (s) for the specific
%   impulse, and (km/s) for the DeltaV.
%
%   EXAMPLE:
%
%       DeltaV = tsjiolkovsky_high(-1, 300, 1000, 150)
%
%           DeltaV =
%
%                    5.58131
%
%       M0 = tsjiolkovsky_high([1.2; 2.4], [175 210], -1, [50 100])
%
%           M0 =
%
%                   100.6105
%                   320.7173
% See also .

% Author:
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 9/Apr/2009.

    % constants
    g0   = 9.80665/1000;  % [km s-2]
    ceff = Isp*g0;        % [km/s]
    
    % calculate output
    if (M0 == -1)    , out = Me(:).*exp(DeltaV(:)./ceff(:));  end
    if (Me == -1)    , out = M0(:)./exp(DeltaV(:)./ceff(:));  end
    if (Isp == -1)   , out = DeltaV(:)./log(M0(:)./Me(:))/g0; end
    if (DeltaV == -1), out = ceff(:).*log(M0(:)./Me(:));      end
end
