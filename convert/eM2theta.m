function th = eM2theta(e, M, tol)  
% EM2THETA      Convert eccentricity and mean anomaly to true anomaly
%
% [th] = EM2THETA(e, M) converts the the eccentricity [e] and mean 
% anomaly [M] to the corresponding true anomaly [th]. The mean anomaly 
% [M] should be given in radians. This conversion is essentially an
% black-box function that solves Kepler's equation for any arbitrary 
% conic section; it works correctly for all types of conic sections, 
% and handles scalar/vector/matrix input.
%
% EM2THETA(e, M, tol) (with a third argument) uses a maximum error 
% tolerance for the iterative procedure given by [tol]. The default 
% is 1e-12.
%
% Algorithm: To calculate [theta] from [M], EM2THETA uses a carefully 
% selected first approximate root of Kepler's equation, followed by a 
% Newton-Raphson iteration scheme if the first estimate is not within 
% the limits set by [tol]. See [Seppo Mikkola, "A cubic approximation 
% for Kepler's equation", 1987] for more details.
%
% See also eM2E, eE2M, eE2theta.


% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology

% Last edited 01/Nov/2009
        
    % set default tolerance        
    if (nargin == 2) || isempty(tol), tol = 1e-12; end

    % initialize       
    par = (e == 1);     % parabolic cases        
    th  = NaN(size(e)); % initialize output argument

    % get eccentric anomalies and calculate theta for 
    % non-parabolic cases
    if ~all(par(:)) 
        Enp = eM2E(e(~par), M(~par), tol);
        th(~par) = eE2theta(e(~par), Enp);        
    end

    % handle parabolic cases
    % (analytic solution to Barker's equation)
    if any(par(:))            
        % extract proper [M]'s
        Mp = M(par);     
        % prevent 1-over-zero warnings
        Mp(Mp == 0) = realmin;
        % compute [th]
        y = abs(atan(1/3./Mp));
        x = atan(tan(y/2).^(1/3));
        th(par) = 2*atan(2./tan(2*x));            
        % make sure special cases are exact
        th(par & M == 0) = 0;
        % make sure the quadrant is also correct
        th(par) = th(par).*sign(Mp);            
    end
        
end
