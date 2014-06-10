% method of collocation
function varargout = collocation(X, type, varargin)

    %% initialize
    
    
    %% optimize
    
    
   
    
    %% objective and constraint functions
    
    % objective function, and its gradient
    function [F, gradF] = objective_function(X)
        F = -X(end-3);
        gradF = zeros(size(X));
        gradF(end-3) = -1;
    end

    % nonlinear constraints
    function [c, ceq, grad_c, grad_ceq] = nonlinear_constraints(X)

        % constants
         muS = params.constants.mu_central;
         g0  = params.constants.g0;
         Isp = params.Isp;

        % extract time of flight
        tof = X(1);   X(1) = [];

        % number of points 
        N = length(X)/10;

        % split matrix in control and state parts
        XX = reshape(X, 10, length(X)/10);
        Xs = XX(1:7, :);    Xc = XX(8:10, :);

        % time intervals
        dt = tof/N;

        % equations of motion
        f = EquationOfMotion(Xs, Xc);

        % Compute values of the center point of each node for first leg
        Xchat  = (Xc(:, 1:end-1) + Xc(:, 2:end))/2;
        Xshat  = (Xs(:, 1:end-1) + Xs(:, 2:end))/2 - dt*diff(f(:, 1:end), 1, 2)/8;
        VXshat = 3*diff(Xs(:, 1:end), 1, 2)/2/dt - (f(:, 1:end-1) + f(:, 2:end))/4;

        % equations of motion in the center of the intervals
        fhat = EquationOfMotion(Xshat, Xchat);

        % inequality constraints
        c = sqrt(sum(Xc.^2)) - maxT./Xs(7, :); % force may not exceed maximum possible
        c = [c, fhat(7, :), VXshat(7, :)];     % all massflows must be negative

        % equality constraints
        ceq = VXshat(:) - fhat(:);             % interpolation must equal the state equation

        % gradients come later
        grad_ceq = [];  grad_c = [];

        % equations of motion
        function deriv = EquationOfMotion(Xs, Xc)

            % initialize derivative
            deriv = zeros(7, size(Xs, 2));

            % assign velocities
            deriv(1:3, :) = Xs(4:6, :);

            % compute accelerations
            r32  = (sum(Xs(1:3, :).^2)).^(3/2);
            rr32 = [r32; r32; r32];
            deriv(4:6, :) = -muS.*Xs(1:3, :)./rr32 + Xc;

            % massflow
            deriv(7, :) = -Xs(7, :).*sqrt(sum(Xc.^2)) /g0/Isp;

        end
    end

end
    