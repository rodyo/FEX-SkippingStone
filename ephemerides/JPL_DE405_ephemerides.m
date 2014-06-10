function statevec = JPL_DE405_ephemerides(times, ppobjects)
% JPL_DE405_EPHEMERIDES         Calculates JPL/DE405 planetary ephemerides
%
%
% See also two_body_ephemerides, quintic_ephemerides, 
% initialize_ephemerides_generators.


% NOTE: the std.grav.param. of the SUN is copy-pasted in the 
% failsafe-cases; I was just too lazy to adjust every function 
% that depended on this one at the time ^_^ In the future, this 
% should of course be generalized to CONSTANTS.mu_central

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology
%              Partly performed at ISAS/JAXA

    % Last Edited: 22/Oct/2009
            
    % all break points are the same
    breaks = ppobjects{1}.breaks; 
     
    % JPL/DE405 ephemerides beyond 2099 are not available for
    % some bodies. In those (rare) cases, generate the ephemerides
    % with PROGRESS_ORBIT()
    if times > breaks(end)
        % first get the last ephemerides in the dataset
        last_statevec = JPL_DE405_ephemerides(breaks(end), ppobjects(:));
        % get the time step required
        timestep = times - breaks(end);
        % generate ephemerides
        statevec = progress_orbit(timestep, last_statevec, 132712439940);
        % we're done
        return
    end
    
    % JPL/DE405 ephemerides before 1/Jan/2000 are not available. Use
    % PROGRESS_ORBIT() for those cases too
    if times < 0
        % first get the first ephemerides in the dataset
        first_statevec = JPL_DE405_ephemerides(0, ppobjects(:));
        % generate ephemerides
        % (time step required is the time itself)
        statevec = progress_orbit(times, first_statevec, 132712439940);
        % we're done
        return
    end
    
    % Everything below was copy-pasted from PPVAL() and adjusted for speed
    
    % evaluate single site. We know beforehand *all* PP's have
    % breakpoints at intervals of 10 days. So we can make use of that:
    ind = max(1, ceil(times/10)); % also make sure the MINIMUM index is [1], not [0]
    % the more general way (much slower for this purpose)
    % ind = nnz([-inf, breaks(2:end-1), inf] <= times)
    
    % shift times to proper evaluation sites (all breakpoints are the same)
    t2 = times - breaks(ind);
    % create time-cubed array
    tt = [t2.^3; t2.^2; t2; 1];
    
    % calculate coordinates and speeds
    x =  ppobjects{1}.coefs(ind, :) * tt;     xdot = ppobjects{4}.coefs(ind, :) * tt;
    y =  ppobjects{2}.coefs(ind, :) * tt;     ydot = ppobjects{5}.coefs(ind, :) * tt;
    z =  ppobjects{3}.coefs(ind, :) * tt;     zdot = ppobjects{6}.coefs(ind, :) * tt;
        
    % insert results statevecs vector
    statevec = [x, y, z, xdot, ydot, zdot];
        
end % JPL DE/405 ephemerides calculator
