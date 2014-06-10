% FINDREALROOTS     Find approximations to all real roots of any function 
%                   on an interval [a, b].
%
% USAGE:
%  Roots = FindRealRoots(funfcn, a, b, n, vectorized, make_plot)
%
% FINDREALROOTS() approximates all the real roots of the function 'funfcn' 
% in the interval [a,b]. It does so by finding the roots of an [n]-th degree 
% Chebyshev polynomial approximation, via the eignevalues of the associated 
% companion matrix. 
%
% When the argument [vectorized] is [true], FINDREALROOTS() will evaluate 
% the function 'funfcn' at all [n] required points in the interval 
% simultaneously. Otherwise, it will use ARRAFUN() to calculate the [n] 
% function values one-by-one. [vectorized] defaults to [false]. 
%
% When the argument [make_plot] is true, FINDREALROOTS() plots the 
% original function and the Chebyshev approximation, and shows any roots on
% the given interval. Also [make_plot] defaults to [false]. 
%
% All [Roots] (if any) will be sorted.
% 
% First version 26th May 2007 by Stephen Morris, 
% Nightingale-EOS Ltd., St. Asaph, Wales.
%
% Modified 14/Nov by Rody Oldenhuis
% Delft university of Technology, Delft.
%
% See also roots, eig.

function Roots = FindRealRoots(funfcn, a, b, n, vectorized, make_plot)

    % parse input and initialize.
    inarg = nargin; 
    if n <= 2, n = 3; end                    % Minimum [n] is 3:    
    if (inarg < 5), vectorized = false; end  % default: function isn't vectorized
    if (inarg < 6), make_plot = false; end   % default: don't make plot
    
    % some convenient variables
    bma = (b-a)/2;  bpa = (b+a)/2;  Roots = [];

    % Obtain the Chebyshev coefficients for the function
    %
    % Based on the routine given in Numerical Recipes (3rd) section 5.8;
    % calculates the Chebyshev coefficients necessary to approximate some
    % function over the interval [a,b]
    
    % initialize 
    c = zeros(1,n);  k=(1:n)';  y = cos(pi*((1:n)-1/2)/n); 
    % evaluate function on Chebychev nodes
    if vectorized
        f = feval(funfcn,(y*bma)+bpa);
    else
        f = arrayfun(@(x) feval(funfcn,x),(y*bma)+bpa);
    end
    
    % compute the coefficients
    for j=1:n, c(j)=(f(:).'*(cos((pi*(j-1))*((k-0.5)/n))))*(2-(j==1))/n; end       
        
    % coefficients may be [NaN] if [inf]
    % ??? TODO - it is of course possible for c(n) to be zero...
    if any(~isfinite(c(:))) || (c(n) == 0), return; end
        
    % Define [A] as the Frobenius-Chebyshev companion matrix. This is based
    % on the form given by J.P. Boyd, Appl. Num. Math. 56 pp.1077-1091 (2006).
    one = ones(n-3,1);
    A = diag([one/2; 0],-1) + diag([1; one/2],+1);
    A(end, :) = -c(1:n-1)/2/c(n);
    A(end,end-1) = A(end,end-1) + 0.5;
    
    % Now we have the companion matrix, we can find its eigenvalues using the
    % MATLAB built-in function. We're only interested in the real elements of
    % the matrix:
    eigvals = eig(A);  realvals = eigvals(imag(eigvals)==0);
    
    % if there aren't any real roots, return
    if isempty(realvals), return; end

    % Of course these are the roots scaled to the canonical interval [-1,1]. We
    % need to map them back onto the interval [a, b]; we widen the interval just
    % a tiny bit to make sure that we don't miss any that are right on the 
    % boundaries.
    rangevals = nonzeros(realvals(abs(realvals) <= 1+1e-5));

    % also sort the roots
    Roots = sort(rangevals*bma + bpa);

    % As a sanity check we'll plot out the original function and its Chebyshev
    % approximation: if they don't match then we know to call the routine again
    % with a larger 'n'.
    if make_plot
        % simple grid
        grid = linspace(a,b, max(25,n));
        % evaluate function
        if vectorized
            fungrid = feval(funfcn, grid);
        else
            fungrid = arrayfun(@(x) feval(funfcn,x), grid);
        end        
        % corresponding Chebychev-grid (more complicated but essentially the same)
        y = (2.*grid-a-b)./(b-a); d = zeros(1,length(grid)); dd = d;
        for j = length(c):-1:2, sv=d; d=(2*y.*d)-dd+c(j); dd=sv; end, chebgrid=(y.*d)-dd+c(1);
        % Now make plot
        figure(1), clf,  hold on
        plot(grid, fungrid ,'color' , 'r');
        line(grid, chebgrid,'color' , 'b'); 
        line(grid, zeros(1,length(grid)), 'linestyle','--')
        legend('function', 'interpolation')
    end % make plot
        
end % FindRealRoots
