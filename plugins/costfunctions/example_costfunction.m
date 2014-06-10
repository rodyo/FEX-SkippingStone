function varargout = example_costfunction(varargin)
% NOTE: don't use "global MainWin" here - this would give difficulties when
% this function is evaluated on a parallel computer cluster

    %% Return info only

    costfun.name = 'Maximum amount of MP''s';
    costfun.description = [...
        'This is en example costfunction that only serves to show ',...
        'how to write custom costfunctions.'];
    costfun.function_handle = @(result) costfcn(result);
    costfun.axis_label = '';
    
    % return information about this cost function when no input arguments
    % have been provided
    if (nargin == 0), varargout{1} = costfun; return, end
    
    %% Return objective & constraint values        
    % (find the amount & quality of MP's)
    function [cost, constraint, output_data] = costfcn(varargin)
        % the function doesn't actually do anything; return dummy values
        cost        = -inf;
        constraint  = 0;
        output_data = [];        
    end
    
end
