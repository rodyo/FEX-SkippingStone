function varargout = max_MPs(MainWin, result)
% NOTE: don't use "global MainWin" here - this would give difficulties when
% this function is evaluated on a parallel computer cluster

    %% Return info only

    % return information about this cost function when no input arguments
    % have been provided
    if (nargin == 0)
        % info structure
        costfun.name = 'Maximum amount of MP''s';
        costfun.description = [...
            'Maximize both the amount [N] and quality [Q] of minor planets that the spacecraft will ',...
            'come close to. The quality [Q] in this context means a low relative speed [Vrel] and ',...
            'a high value for [model.MPs.scivalue. Note that this cost function is very costly to ',...
            'compute, so it is advised to only use this costfunction if several days of computation ',...
            'time are acceptable, or if the optimization is carried out on a parallel cluster.'];
        costfun.function_handle = @(MainWin,result) max_MPs(MainWin, result);
        costfun.axis_label = 'MP flyby quality [= sum(Qsci/|Vrel|)]';
        % return argument
        varargout{1} = costfun; return
    end
    
    %% Return objective & constraint values        
    % (find the amount & quality of MP's)
    
    % this costfunction does not have any constraints associated with it
    constraints = [];    
    % initially, there are also no additional data
    output_data = [];
    
    % when patched-conics result has constraint violations
    if result.is_violated
        cost = 0;
       
    % when patched-conics result is fully feasible
    else  
        
        % NOW get relevant data
        model     = getappdata(MainWin, 'model'    );
        settings  = getappdata(MainWin, 'settings' );
        constants = getappdata(MainWin, 'constants');
        
        % high thrust or low thrust?
        high_thrust = settings.propulsion.selected(1);
        low_thrust  = ~high_thrust;
        
        % include scientific value in the result?
        use_scivalue = isfield(model.MPs, 'scivalue'); 
        
%??? TODO - generalize
mindist_options = struct(...
    'threshold', [0.01 0.09 30]*constants.AU);
        
        % convert transfer and initial times to something more convenient
        times = cumsum([result.t0, result.tfs]);
        
        % initialize loop
        reachable = false(size(model.MPs.as,1),numel(result.tfs));% logical        
        encounter_times = zeros(size(reachable));% double
        min_dist        = encounter_times;       % double
        rel_speed       = encounter_times;       % double
        if use_scivalue
            scivalue = repmat(model.MPs.scivalue, 1, numel(result.tfs)); % double
        end
                 
        % loop through all trajectories in the result
        for jj = 1:numel(result.tfs)
                        
            % put all MP's at the proper initial positions
            % NOTE: don't forget to convert times to seconds!
            t0 = times(jj); tend = times(jj+1);
            MPs_new_Ms = model.MPs.M0s + ...
                model.MPs.ns.*(t0-model.MPs.epoch_days_past_J2000)*86400;
            MPs_new_theta0s = eM2theta(model.MPs.es, MPs_new_Ms);
            
            % run the minimum distance routine
            % (high-thrust post-processor differs from low-thrust post-processor)
            if high_thrust
                % convert all trajectories therein to orbital elements
                sc_elements = cart2kep([result.r_departure, result.V_departure], ...
                    constants.mu_central, 'theta');
                % run routine  
                [reachable(:,jj), encounter_times(:,jj), min_dist(:,jj), rel_speed(:,jj)] = ...
                    minimum_distance_conics(sc_elements(jj,:), ...
                    [model.MPs.as, model.MPs.es, model.MPs.is, ...
                    model.MPs.Omegas model.MPs.omegas, MPs_new_theta0s], ...
                    t0, tend, constants.mu_central, mindist_options);                
            elseif low_thrust
                % run routine
                [reachable(:,jj), encounter_times(:,jj), min_dist(:,jj), rel_speed(:,jj)] = ...
                    minimum_distance_exposins(result.exposins(jj,:), result.R{jj},...
                    [model.MPs.as, model.MPs.es, model.MPs.is, ...
                    model.MPs.Omegas, model.MPs.omegas, MPs_new_theta0s], ...
                    t0, tend, constants.mu_central, mindist_options);
            end
            
        end % minimum-distance for loop
                     
        % gather data
        reachables = any(reachable,2);
        output_data.reachable  = reachable;
        output_data.min_dist   = min_dist(reachables,:);
        output_data.rel_speed  = rel_speed(reachables,:);
        output_data.encounter_times = encounter_times(reachables,:);
        output_data.nearby_MPs = nnz(reachables);
        if use_scivalue, output_data.scivalue = scivalue(reachables,:); end
        
        % compute costfunction (including scientific value)
        if use_scivalue
            % compute quality measure
            Q = output_data.scivalue./output_data.rel_speed;
%??? TODO - If specific MP's are found more than once, this
% method will fail
Q(~isfinite(Q)) = 0; Q = sum(Q(:));
            % value of the cost function
            cost = -Q;
            
            % compute costfunction (excluding scientific value)
        else
            cost = -output_data.nearby_MPs;
        end
                
    end % if violated (or not)
    
    % assign output arguments
    varargout{1} = cost;
    varargout{2} = constraints;
    varargout{3} = output_data;
    
end
