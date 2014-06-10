% NBODY        N-body integrator for [M] insignificant bodies

function [t, states, output] = NBody_Battin(y0, tspan, options)

    %% initialize
        
    % basic error trap
    error(nargchk(3,3,nargin));
    
    % extract all data (the fool-proof way)
    
    % ephemerides generator
    if ~isfield(options, 'ephemerides') || ~isa(options.ephemerides, 'function_handle')
        errordlg(['OPTIONS.ephemerides should contain a function-handle to \n',...
            'the ephemerides generator for the perturbing bodies.'], 'ERROR');
    else
        ephem_generator = options.ephemerides;
    end
    % standard gravitational parameters
    if ~isfield(options, 'GMs')
        
    else
        GMs = options.GMs;
    end
    % integrator
    if ~isfield(options, 'integrator')
        
    else
        integrator = options.integrator;
    end
    % perturbing bodies
    if ~isfield(options, 'perturbers')
        
    else
        N = options.perturbers;
    end    
    
    % even more error traps
    if isempty(N) || (isscalar(N) && N == 1)
        warning('NBody:problem_is_Keplerian', ...
            ['The given problem is equivalent to a two-body problem; \n',...
            'Results can be obtained more quickly by using Kepler''s equation.']);        
    end
    if numel(GMs)>numel(N)
        warning('Nbody:N_smaller_GMs',...
            ['List of gravitational parameters is longer than the list\n',...
            'of perturbing bodies. Is this correct?'])
    end
    if numel(GMs) < numel(N)
        error('Nbody:N_larger_GMs',...
            ['More perturbing bodies were given than gravitational parameters.\n',...
            'NBODY() cannot continue.'])        
    end
    if ~any(size(y0) == 6)
        error('NBody:inconsistent_dimensions',...
            'The initial statevectors [y0] should have length 6 in at least one dimension.')
    end
    
    % make sure the central body (body #1) is in N
    N = sort(N(:).'); if ~any(N==1), N = [1,N]; end        
    % also make sure [y0] has the right size (columns are more efficient)
    if size(y0,2)==6, y0 = y0.'; end        
    % amount of insignificant bodies being integrated
    M = size(y0,2);
    
tspan = tspan * 86400;
    
    % initialize all fields in output
    output.ephemerides     = 0;
    output.ODEevaluations  = 0;
    output.sucessful_steps = 0;
    
    %% CALL INTEGRATOR
    
    % NOTE:
    % the velocities are smallest in magnitude; a tolerance of 1e-3 means a
    % maximum error in the velocity of 1 m/s/step, which corresponds to a 
    % position error of 1 m/step. This is a bit superfluous and unneccesarily 
    % slow; a better approach would be to use some sort of scaling/re-scaling 
    % to make the integration run faster, without loosing accuracy in both
    % velocity and position.
    options = odeset('abstol', 1e-3); 
    
    % Integrators may be of type Runge/Kutta (ODE); these can run as-is.
    try 
        ode = true; [t, y] = integrator(...
            @(t,y) dydt(t,y,ode), tspan, y0(:), options);
            
    % Integrators may also be of type Runge/Kutta/Nystrom (RKN). In
    % this case, an additional initial estimate (dydt(t0)) is needed.
    catch%#ok   
        try % this may also fail of course...
            % first split the given input
            dydt0 = y0(4:6,:);    y0 = y0(1:3,:);
            % integrate (with the formula for the *second* derivative)
            ode = false; [t, y, dy] = integrator(...
                @(t,y) dydt(t,y,ode), tspan, y0(:), dydt0(:), options);
        catch ME
            ME2 = MException('NBody:integrator_failure',...
                'Integrator failed to evaluate. NBody cannot continue.');            
            rethrow(addCause(ME,ME2));
        end
    end
    
    %% Parse output
    
    % [y] (and [dy]) are concatenated, so split the bodies up 
    % again and insert each one in a more convenient cell-array
    states = cell(M,1);
    for i = 1:M
        if ode
            states{i} = y(:,1:6);  y = y(:, 7:end);
        else            
            states{i} = [y(:, 1:3), dy(:, 1:3)];  
            y = y(:, 4:end); dy = dy(:, 4:end);
        end
    end  
    
    % number of sucessful steps 
    output.sucessful_steps = size(states{1},1);
        
    %% the derivative
    
    % Both can be done using Battin's reformulation, which reduces roundoff
    % error and thus increases the accuracy per step
    
    % computation of [y'] (for ODE-type methods) or [y''] (for RKN-type methods)
    function F = dydt(t, y, ode)
                 
        % first reshape input
        y = reshape(y,3+3*ode,[]);  
        % initialize output
        F = zeros(size(y));
        
        % absolute location of the central body at this time
        stateC    = ephem_generator(1, t/86400); 
        posC(:,1) = stateC(1:3);
        
        % keep track of ephemerides calculations
        output.ephemerides = output.ephemerides + 1;
        
        % first compute accelerations of insignificant bodies w.r.t. central body
        r   = bsxfun(@minus, y(1:3,:), posC); % vectors from insignificant bodies to central body
        rm3 = sum(r.*r).^(1.5);               % their magnitude, to the third power        
        % accelerations due to central body alone
        acc = -GMs(1)*bsxfun(@rdivide, r, rm3);
        
        % loop through all the perturbing bodies
        for j = N(2:end)
            % calculate the positions of this perturbing body at this time
            states = ephem_generator(j, t/86400); posj(:,1) = states(1:3); 
            % keep track of ephemerides calculations
            output.ephemerides = output.ephemerides + 1;
            % compute difference vectors
            rhoj = posj - posC;             % central body to perturbing body
            dj   = bsxfun(@minus, r, posj); % perturbing body to insignificant bodies            
            % compute their magnitudes, to the correct power            
            rhojm = sqrt(rhoj.'*rhoj);
            djm3  = (sum(dj.*dj)).^(1.5) ;           
            % compute f(q)
            qj  = sum( (r/rhojm/rhojm) .* bsxfun(@minus, r, 2*rhoj) );
            fqj = qj .* ( (3+qj.*(3+qj)) ./ (1+(1+qj).^(3/2)) );        
            % Battin's re-formulation (MATLAB'ed for an amount of M-small bodies)          
            acc = acc - GMs(j)*bsxfun(@rdivide, r + bsxfun(@times, fqj, rhoj), djm3);
        end
        
        % complete the derivative
        F((1:3)+3*ode,:) = acc; 
        if ode, F(1:3,:) = y(4:6,:); end 
        F = F(:);   
        
        % keep track of the amount of derivative-evaluations
        output.ODEevaluations = output.ODEevaluations + 1;
        
    end % yprime
    
end % N-body integrator for M insignificant bodies
