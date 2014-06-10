% PATCHED_CONICS
%
%
%
%
function result = patched_conics(seq, X, params)
    
    %% Initialize
    
    % the format of [X] depends on the solution type
    % =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    % FORMAT (high-thrust/exposins/equinoctial elements):
    %
    %     [t0,...
    %     tf(1)  longway(1)  k2(1)   m(1)   leftbranch(1),...
    %     tf(2)  longway(2)  k2(2)   m(2)   leftbranch(2),...
    %     ...
    %     tf(N)  longway(N)  k2(N)   m(N)   leftbranch(N)]
    %
    % with [tf] > 0 the times of flight, [longway] = +- 1 the
    % corresponding long-way solution, [m] >= 1 the amount of 
    % complete revolutions to use, [k2] > 0.01 the free 
    % Lambert parameter for ExpoSin trajectories, and 
    % [leftbranch] = +-1 the selector for using the left-or 
    % right branch (in case [m] > 0)
    % =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    % =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    % FORMAT (collocation/Sins & Flanagan/Ballistic integration):
    %
    %     [t0,...
    %     tf(1), tf(2), ..., tf(N),...
    %     <other data>]
    %
    % 
    % =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    % launch date is always the first element in [X]
    t0 = X(1);
    
    % the format of [X] depends on the solution type:
    
% X(2) =[];
    
    % Second-order methods        
    if any(strcmpi(params.solution_type, {'Sims & Flanagan';'collocation';'ballistic-integration'}))        
        % the first [number_of_targets]-elements are the times of flight        
        tfs = X(2:numel(seq)); tfs = tfs(:).';
        % what the other entries represent depends on the type of method
        Y = X((numel(seq)+1):end);
        % the other data really depends on the specific method
        switch lower(params.solution_type)
            case 'sims & flanagan'
                N = params.N; % number of steps
                
            case 'collocation'
                N = params.N; % number of points
                
            case 'ballistic-integration'
                % CBF-sequences: extract location of pericenter and Delta-V, in 
                % addition to departure velocities at all GAM-bodies
                CBF_indices = strcmpi(params.types, 'Central Body Flyby');                                
                if any(CBF_indices)
%??? TODO: GENERALIZE                     
% pericenter location and DeltaV
DV_CBF = Y(1);
% departure velocities at GAM-bodies
V_init = reshape(Y(2:end),3,[]).';
%??? TODO: GENERALIZE 
                                       
                % non-CBF sequences: just departure velocities at GAM-bodies                
                else
                    V_init = reshape(Y,3,[]).';
                end                
        end
        
        % use dummies for everything else
        longway = [];   k2s = [];   ms = [];   leftbranch = [];
            
    % First-order methods
    else
        % extract TOF, longway, m, and DSM parameters from X
        tfs        = X(2:5:end);        % times of flight
        longway    = sign(X(3:5:end));  % use long- or short way?
        k2s        = X(4:5:end);        % ExpoSin - parameters
        ms         = round(X(5:5:end)); % be careful to make integers of [ms]
        leftbranch = sign(X(6:5:end));  % take left or right branch?
    end
    
    
    % some initializations
    number_of_targets = numel(tfs);  % total number of bodies
    ii = 1;                          % changeable loop-index 
    
    % rename solvers 
    ephemerides = params.ephemerides_generator;       
    
    % make absolutely sure there can be *NO* Delta-V 
    % for completely un-powered trajectories
    if all(strcmpi(params.types, 'un-powered'))
        params.max_total_DV = 0;    params.max_GAM_DeltaV(:) = 0;
    end
    
    % declare some zero/NaN arrays
    NaN_1  = NaN(number_of_targets-1, 1);     NaN_2  = NaN(number_of_targets, 3);
    zero_1 = zeros(number_of_targets-1, 1);   zero_2 = zeros(number_of_targets, 1);
    
    % declare output to be calculated      
    result.seq              = seq;          result.r_departure         = NaN_2;
    result.rps              = NaN_1;        result.r_arrival           = NaN_2;
    result.closest_approach = NaN_1;        result.V_departure         = NaN_2;
    result.Vinf_launch      = NaN(1, 3);    result.V_arrival           = NaN_2;
    result.Vinf_target      = NaN(1, 3);    result.V_departure_GAMbody = NaN_2;
    result.Vinf_plus        = NaN_2;        result.V_target_GAMbody    = NaN_2;
    result.Vinf_min         = NaN_2;        result.turnangle           = NaN_1;
    result.DeltaVs          = NaN_1;        result.additional_DeltaV   = NaN_1;    
    result.DeltaVtot        = 0;            result.masses              = [NaN_1;NaN;NaN];
    result.endmass          = NaN;          result.longway             = longway;
    result.m                = ms;           result.tfs                 = tfs; 
    result.t0               = t0;           result.leftbranch          = leftbranch;
    result.k2_parameters    = k2s;          result.mindist_to_central  = NaN_1;
    result.maxdist_to_central = NaN_1;      result.C3                  = inf;
    
    % declare constraint violations
    V.C3_violation          = 0;            V.rp_violation       = zero_1;
    V.GAM_DV_violation      = zero_1;       V.mindist_violation  = zero_2;
    V.TOF_violation         = 0;            V.GAM_Vinf_violation = zero_1;
    V.total_DV_violation    = 0;            V.endmass_violation  = 0;
    V.Vinf_target_violation = 0;            V.CBF_violation      = 0;
    V.miss_distance         = zero_2;       V.other_constraints  = [];    
    result.violations       = V;            result.is_violated   = false;
    
    % results specific to the ExpoSin-lambert targeter 
    if strcmpi(params.solution_type, 'exposins')
        % the exposin parameters
        result.exposins = NaN(number_of_targets, 10);
        % associated rotation matrix
        result.R        = cell(number_of_targets,1);
    end
    
    % results specific to Sims & Flanagan/collocation methods
    if any(strcmpi(params.solution_type, {'Sims & Flanagan';'collocation'}))
        % number of points per leg
        result.N = N;
        % DeltaVs is more conveniently expressed as cell-array
        result.DeltaVs = cell(number_of_targets, 1);
    end
    
    % statistics
    result.ephemerides      = 0;           result.rp.failure           = 0;    
    result.rp.iterations    = 0;           result.lambert.failure      = 0; 
    result.lambert.success  = 0;           result.lambert.impossible   = 0;
    result.lambert.overflow = 0;           result.central_body.success = 0;          
    result.central_body.impossible = 0;    result.central_body.failure = 0;    
    result.ODEevaluations   = 0;           result.integration_steps    = 0;
    
    % quick exit
    if any(~isfinite(X)), result = all_violated(result); return, end
           
    % cumulative times (for the ephemerides generator)
    times = cumsum([t0, tfs]);
    
    % initial mass
    M = params.M0; result.masses(1) = M;
    
    % initially, NO central body flybys
    cbf = [];
    
    % determine ephermeris for the first body
    xs = ephemerides(seq(1), times(1));    
    % keep track of the amount of ephemerides calculations
    result.ephemerides = 1;    
    % extract vectors from ephemerides
    r_departure         = xs(1, 1:3); % departure coordinates and velocity     
    V_departure_GAMbody = xs(1, 4:6); % of the departure body
             
    % save these vectors          
    result.r_departure(1, :) = r_departure;
    result.V_departure_GAMbody(1, :) =  V_departure_GAMbody;
                    
    %% LOOP THROUGH THE SWINGBYS
            
    % Doing swingby's one-by-one is a bit slower but a LOT easier to program
    while ii <= number_of_targets
            
        %% SKIP CENTRAL BODY FLYBYS 
        
        if (ii < number_of_targets) && strcmpi(params.types(ii), 'Central Body Flyby')
            
            % first-order ballistic 
            if strcmpi(params.solution_type, 'ballistic')
                % store this [ii]
                cbf = [cbf; ii];%#ok
                % get ephemerides of the terminal GAM-body
                xs = ephemerides(seq(ii+2), times(ii+2));
                result.ephemerides = result.ephemerides + 1;
                % extract vectors from ephemerides
                r_target = xs(1:3);  V_target_GAMbody = xs(4:6);
                % insert arrival location in result
                result.r_arrival (ii+1, :) = r_target;
                result.V_target_GAMbody(ii+1, :) = V_target_GAMbody;
                if (ii+2 <= number_of_targets)
                    result.r_departure (ii+2, :) = r_target;
                    result.V_departure_GAMbody(ii+2, :) = V_target_GAMbody;
                end
                % copy quantities for next iteration
                r_departure = r_target;
                V_departure_GAMbody = V_target_GAMbody;
                % increase loop index
                ii = ii + 2;
                % and continue
                continue;
                
            % second-order ballistic
            elseif strcmpi(params.solution_type, 'ballistic-integration')

% get ephemerides of the terminal GAM-body
xs = ephemerides(seq(ii+2), times(ii+2));                
% extract vectors from ephemerides
r_target = xs(1:3);  V_target_GAMbody = xs(4:6);
% Save next planetary velocity, and next target radius   
if (ii < number_of_targets)
    result.V_departure_GAMbody(ii+1, :) = V_target_GAMbody;
    result.r_departure (ii+1, :) = r_target;
end

% Copy
V_departure = V_init(ii,:);
% Lambert problem ( [state0], [tspan] )
[~, states] = params.integrator(...
    [r_departure +  1e4*V_departure, V_departure], times(ii)+[0 tfs(ii)]);

% copy quantities for next iteration
r_departure = r_target;
V_departure_GAMbody = V_target_GAMbody;

% insert  positions & velocities in result
final_state = states{1}(end,:);
r_CBF  = final_state(1:3);
V_CBF  = final_state(4:6);

% save results
result.r_arrival(ii,:) = r_CBF;
result.r_departure(ii+1,:)= r_CBF;
result.V_arrival(ii,:) = V_CBF;
% add Delta-V
V_CBF_unit = V_CBF / sqrt(V_CBF*V_CBF.');
V_CBF = V_CBF + DV_CBF*V_CBF_unit;
result.DeltaVs(ii) = DV_CBF;
result.additional_DeltaV(ii) = 0;
% save results
result.V_departure(ii+1,:) = V_CBF;

% Lambert problem ( [state0], [tspan] )
[~, states] = params.integrator(...
    [r_CBF, V_CBF], times(ii+1)+[0 tfs(ii+1)]);

% insert  positions & velocities in result
final_state = states{1}(end,:);
r_terminal  = final_state(1:3);
V_arrival   = final_state(4:6);

result.V_arrival(ii+1,:) = V_arrival;

% calculate the miss distance
miss_vector = (r_terminal - r_target);
result.miss_distance(ii+1) = sqrt(miss_vector*miss_vector.');

% increase loop index
ii = ii + 2;
% and continue
continue;


                
            end
        end
        
        %% GET EPHEMERIDES AND SOLVE LAMBERT PROBLEM
            
        % determine ephermeris for the next target body
        xs = ephemerides(seq(ii+1), times(ii+1));
        result.ephemerides = result.ephemerides + 1;
        % extract vectors from ephemerides
        r_target = xs(1:3);  V_target_GAMbody = xs(4:6);
        % insert arrival location in result
        result.r_arrival (ii, :) = r_target;
        result.V_target_GAMbody(ii, :) = V_target_GAMbody;
        % Save next planetary velocity, and next target radius   
        if (ii < number_of_targets)
            result.V_departure_GAMbody(ii+1, :) = V_target_GAMbody;
            result.r_departure (ii+1, :) = r_target;
        end
         
        % Call Lambert-targeter
        switch lower(params.solution_type)
            % hight-thrust ballistic flight, two-body approximation. The Lambert-targeter 
            % returns the terminal velocities and the minimum and maximum distances to 
            % the central body.
            case 'ballistic'
                % Lambert problem
                [V_departure, V_arrival, extremal_distances, exitflag] = lambert_high(...
                    r_departure, r_target, longway(ii)*tfs(ii), leftbranch(ii)*ms(ii), params.muC);
                % insert in result
                result.mindist_to_central(ii) = extremal_distances(1);
                result.maxdist_to_central(ii) = extremal_distances(2);
                % rename the minimum distance
                mindist_to_central = extremal_distances(1);
                                
            % hight-thrust ballistic flight, numerical integration. The "Lambert-targeter"
            % returns the complete statevector at each integration step, but only the 
            % terminal one is needed here:            
            case 'ballistic-integration'
                % Copy                 
                V_departure = V_init(ii,:);
                r_departure = r_departure + 1e4*V_departure; % ??? USE OFFSET
                % Lambert problem ( [state0], [tspan] )
                [~, states, output] = params.integrator(...
                    [r_departure, V_departure], times(ii)+[0 tfs(ii)]);
                
%states
result.DETAIL_LEG{ii} = states;
                
                % insert  positions & velocities in result
                final_state = states{1}(end,:);                
                V_arrival   = final_state(4:6);
                r_terminal  = final_state(1:3);
                % calculate the miss distance
                miss_vector = (r_terminal - r_target);
                result.miss_distance(ii) = sqrt(miss_vector*miss_vector.');                
                % the minimum distance to the central body can't be computed easily; so 
                % just set it equal to the minimum allowable distance so that this
                % constraint is never violated
                mindist_to_central = params.mindist_to_central;
                % keep track of the statistics    
                result.ODEevaluations    = result.ODEevaluations + output.ODEevaluations;
                result.ephemerides       = result.ephemerides + output.ephemerides;
                result.integration_steps = result.integration_steps + output.sucessful_steps;
                %??? TODO
                exitflag = 1;
                        
            % Exponential Sinusoids
            % The lambert targeter returns the terminal velocities, the minimum and maximum 
            % distance to the central body, and the estaimted mass left after flying the 
            % ExpoSin-trajectory.
            case 'exposins'                
                % Lambert problem
                [V_departure, V_arrival, extremal_distances, M, ...
                    result.exposins(ii,:), result.R{ii}, exitflag] = lambert_low_exposins(...
                    r_departure, r_target, longway(ii)*tfs(ii), k2s(ii), leftbranch(ii)*ms(ii), M, ...
                    params.Isp, params.muC);
                % save new endmass 
                result.masses(ii+1) = M;                
                % insert in result
                result.mindist_to_central(ii) = extremal_distances(1);
                result.maxdist_to_central(ii) = extremal_distances(2);
                % rename the minimum distance
                mindist_to_central = extremal_distances(1);
                
            case 'equinoctial elements'
                % ??? TODO
                
            case 'sims & flanagan'
                % ??? TODO
                
            case 'collocation'
                % ??? TODO
                
        end
        
        % save results
        result.V_departure(ii, :) = V_departure;
        result.V_arrival  (ii, :) = V_arrival;        
        % keep track of the total amount of lambert problems solved,
        % and save the individual failure modes
        if exitflag == 1 % sucessful
            result.lambert.success = result.lambert.success + 1;
        else % unsucessful
            if exitflag == -1 % no solutions possible
                result.lambert.impossible = 1;
            elseif exitflag == -2 || exitflag == -4  % solutions possible, but algorithm failure
                result.lambert.failure = 1;
            elseif exitflag == -3 % numerical over- or underflow
                result.lambert.overflow = 1;
            end
            % set all constraints to inf and return
            result = all_violated(result); return
        end
        
        % handle minimum distance violation
        result.violations.mindist_violation(ii) = params.mindist_to_central - mindist_to_central;
        
        %% FIRST ITERATION: LAUNCH
        
        % launch
        if (ii == 1)            
            % get launch conditions
            [result.C3, result.Vinf_launch, result.violations.C3_violation] = launch(...
                V_departure, V_departure_GAMbody, ...
                params.max_C3*(1 + params.C3LoverD_tolerance/100));            
            % copy launch V-inf 
            result.Vinf_plus(1,:) = result.Vinf_launch;              
            % save departure vectors for the next iteration        
            r_departure = r_target;   
            V_departure_GAMbody = V_target_GAMbody;  
            % done for this iteration
            ii = ii + 1; continue;
        end    
        
        %% SUBSEQUENT ITERATIONS: SWINGBYS
        
        % extract the previous V_arrival
        V_arrival = result.V_arrival(ii-1, :);
        
        % do the swingby
        [DeltaV, additional_DeltaV, rp, Vinfin, Vinfout, turnangle, rp_violation,...
        GAM_DV_violation, GAM_Vinf_violation, iterations, failure] = swingby(...
        params, V_arrival, V_departure, V_departure_GAMbody, ii);
    
        % statistics
        result.rp.iterations = result.rp.iterations + iterations;

        % swingbys might fail
        if failure, result.rp.failure = 1; result = all_violated(result); return, end
        
        % handle constraint violations
        result.violations.GAM_Vinf_violation(ii-1) = GAM_Vinf_violation;
        result.violations.GAM_DV_violation(ii-1)   = GAM_DV_violation;
        result.violations.rp_violation(ii-1)       = rp_violation;
        
        % save results
        result.DeltaVs          (ii-1, :) = DeltaV;
        result.rps              (ii-1, :) = rp;
        result.closest_approach (ii-1, :) = rp - params.mean_Radii(ii-1);
        result.Vinf_min         (ii-1, :) = Vinfin;
        result.Vinf_plus        (ii  , :) = Vinfout;
        result.additional_DeltaV(ii-1, :) = additional_DeltaV;        
% %??? THESIS STUFF
% if isnan(result.additional_DeltaV(ii-1, :))
% result.additional_DeltaV(ii-1, :) = additional_DeltaV;
% else
%     result.additional_DeltaV(ii-1, :) = result.additional_DeltaV(ii-1, :)+additional_DeltaV;
% end
% %??? THESIS STUFF
        result.turnangle        (ii-1, :) = turnangle;
                
        %% LAST ITERATION: ARRIVAL CONDITIONS
        
        % compute last Vinf and constraint violation
        if (ii == number_of_targets)             
            % ??? TODO - insert max Vinf into params
            [result.Vinf_target, result.violations.Vinf_target_violation] = ...
                target(result.V_arrival(end,:), V_target_GAMbody, 1e5);% just a large value for now
            % copy Vinf-target
            result.Vinf_min(end,:) = result.Vinf_target;
        end 
        
        %% LOOP CONTINUATION
        
        % save radius vectors and mass for the next iteration
        r_departure = r_target;    V_departure_GAMbody = V_target_GAMbody; 
        
        % increase loop index
        ii = ii + 1;
        
    end  % main loop
        
    %% CENTRAL BODY FLYBYS
    
    % are there any CBF-problems to be solved?
    if ~isempty(cbf), result = CBF(params, cbf, result); end
    
    %% FINALIZE SOLUTION
    
    % total DeltaV
    result.DeltaVtot = sum(result.DeltaVs) + sum(result.additional_DeltaV);    
    % total DeltaV might be larger than the maximum possible DeltaV. 
    result.violations.total_DV_violation = result.DeltaVtot - params.max_total_DV;
    
    % total time of flight might be longer than the allowed value. 
    result.violations.TOF_violation = sum(tfs) - params.max_TOF;  
    
    % calculate the mass (for hight-thrust problems)
    if any(strcmpi(params.solution_type, {'ballistic';'ballistic-integration'}))
        % mass after each GAM
        for ii = 1:numel(result.DeltaVs)
            % calculate mass
            M = Tsjiolkovsky(result.DeltaVs(ii)+result.additional_DeltaV(ii), params.Isp, M, -1);
            % insert in result
            result.masses(ii+1) = M;
        end
        % copy end-mass for consistency with low-thrust patched exposins
        % ??? TODO: end-mass follows from Delta-V at target! capture?
        % insertion? entry? -> Costfunction!
        result.masses(end) = result.masses(end-1);
    end
        
    % save the end mass
    result.endmass = result.masses(end);    
    % if the is violated, compute the amount of violation for the final mass 
    result.violations.endmass_violation = params.Me - result.endmass;        
   
    % boolean for violated solutions
    fields = fieldnames(result.violations);
    for ii = 1:numel(fields)
        if any(result.violations.(fields{ii}) > 0), result.is_violated = true; break; end
    end 
                
end

%% Helper functions

% Central body flyby's
function result = CBF(params, cbf, result)
    
    %% loop through all central body flyby problems
    
    for ii = 1:numel(cbf)
        
        % extract current index
        jj = cbf(ii);
        
        % re-calculate number_of_targets
        number_of_targets = numel(result.seq)-1;
        
        % get the positions and velocities for the two GAM bodies, and
        % the arrival/departure velocities when applicable
        r1 = result.r_departure(jj, :);
        r2 = result.r_arrival(jj+1, :);        
        if (jj > 1),                   VLam1 = result.V_arrival(jj-1, :); end
        if (jj+1 < number_of_targets), VLam2 = result.V_departure(jj+2, :); end
        VGAM1 = result.V_departure_GAMbody(jj, :);
        VGAM2 = result.V_target_GAMbody(jj+1, :);        
        
        % select proper functions to use for the Delta-V equation
        
        % launch
        if (jj == 1)
            confcn_r1 =  @(V_departure) ...
                CBF_launch_wrapper(params, V_departure, VGAM1);
            % normal swingby
        else
            confcn_r1 = @(V_departure) ...
                CBF_swingby_wrapper(params, VLam1, V_departure, VGAM1, jj);
        end
        % arrival
        if (jj+1 == number_of_targets)
            confcn_r2 = @(V_arrival) ...
                CBF_target_wrapper(params, V_arrival, VGAM2);
        % normal swingby
        else
            confcn_r2 = @(V_arrival) ...
                CBF_swingby_wrapper(params, V_arrival, VLam2, VGAM2, jj+2);
        end
        
        % solve central body flyby problem
        [V1, V2, output, exitflag] = ...
            central_body_flyby(r1, r2, result.longway(jj)*(result.tfs(jj)+result.tfs(jj+1)), 0, ...
            params.min_alts(jj), confcn_r1, confcn_r2, params.muC, params.scope); 
                
        % insert results
        if exitflag == 1 || exitflag == -2
            % successfully found solution (with possible constraint violations)
            
            % set correct times of flight
            result.tfs(jj) = output.tf1;
            result.tfs(jj+1) = output.tf2;
            
            % handle constraint
            result.violations.GAM_DV_violation(jj) = output.DeltaV - params.max_GAM_DeltaV(jj);
            result.violations.CBF_violation = ...
                max(result.violations.CBF_violation, output.constrviolation);
              
            % insert values
            result.closest_approach(jj) = output.rp - params.mean_Radii(jj);
            
            % statistics
            result.central_body.success = result.central_body.success + 1;
            
            % CBF departure body is patched conics' launch body
            if (jj == 1)
                % DeltaV at central body
                result.DeltaVs(1) = output.DeltaV;
                % set additional DeltaV to zero
                result.additional_DeltaV(1) = 0;
                % evaluate launch() to get proper values
                [result.C3, result.Vinf_launch, result.violations.C3_violation] = launch(V1, VGAM1, ...
                    params.max_C3*(1 + params.C3LoverD_tolerance/100));
                % Save Vinf_departure & arrival
                result.V_arrival(1,:) = output.Vp_vec(1, :); 
                result.V_departure(1,:) = V1;
                % Save Vinf_min and plus
                result.Vinf_plus(1,:) = result.Vinf_launch;
                result.Vinf_min(1,:) = output.Vp_vec(1, :); 
                % insert location of pericenter
                result.r_arrival(1,:)   = output.rp_vec;
                result.r_departure(2,:) = output.rp_vec;
                
            % CBF departure body means normal swingby
            else
                % evaluate swingby() to get proper values
                [result.DeltaVs(jj-1), result.additional_DeltaV(jj-1), result.rps(jj-1), ...
                    result.Vinf_min(jj-1,:), result.Vinf_plus(jj-1,:), result.turnangle(jj-1),...
                    rp_violation, GAM_DV_violation, GAM_Vinf_violation, iterations, failure] = ...
                    swingby(params, VLam1, V1, VGAM1, jj);
                % closest approach (only one not assigned)
                result.closest_approach(jj-1) = result.rps(jj-1) - params.mean_Radii(jj-1);
                
                % statistics
                result.rp.iterations = result.rp.iterations + iterations;
                
                % handle constraint violations
                result.violations.GAM_DV_violation(jj-1) = GAM_DV_violation;
                result.violations.GAM_Vinf_violation(jj-1) = GAM_Vinf_violation;
                result.violations.rp_violation(jj-1) = rp_violation;
                
                % failure
                if failure, result = all_violated(result); return; end
                                
                % set additional DeltaV to zero
                result.additional_DeltaV(jj) = 0;
                
                % set V_arrival, departure
                result.V_arrival(jj,:)   = output.Vp_vec(1, :); 
                result.V_departure(jj,:) = V1;
                
                % set Vinf_min, plus
                result.Vinf_min(jj,:)  =  output.Vp_vec(1, :); 
                result.Vinf_plus(jj,:) =  output.Vp_vec(2, :);
                
                % DeltaV at central body
                result.DeltaVs(jj) = output.DeltaV;                
                % insert location of pericenter
                result.r_arrival(jj,:)   = output.rp_vec;
                result.r_departure(jj+1,:) = output.rp_vec;
            end
            
            % CBF target body is pathced conics target body
            if (jj+1 == number_of_targets)
                % DeltaV at central body                
                result.DeltaVs(end) = output.DeltaV;      
                % set additional DeltaV to zero
                result.additional_DeltaV(end) = 0;
                % evaluate target() to get proper values
                [result.Vinf_target, result.violations.Vinf_target_violation] = ...
                    target(V2, VGAM2, inf);%??? TODO - insert max Vinf here
                % terminal velocities
                result.V_departure(end,:) = output.Vp_vec(2, :); 
                result.V_arrival(end,:)   = V2;  
                % Save Vinf_plus, min
                result.Vinf_plus(end,:) = output.Vp_vec(2, :); 
                result.Vinf_min(end,:) = result.Vinf_target;
                % insert location of pericenter
                result.r_arrival(end-1,:) = output.rp_vec;
                result.r_departure(end,:) = output.rp_vec;
                
            % CBF target body is normal swingby
            else
                % evaluate swingby() to get proper values
                [result.DeltaVs(jj+1), result.additional_DeltaV(jj+1), result.rps(jj+1), ...
                    result.Vinf_min(jj+1,:), result.Vinf_plus(jj+1,:), result.turnangle(jj+1,:),...
                    rp_violation, GAM_DV_violation, GAM_Vinf_violation, iterations, failure] =...
                    swingby(params, VLam2, V2, VGAM2, jj+2);                
                % closest approach (only one not assigned)
                result.closest_approach(jj+1) = result.rps(jj+1) - params.mean_Radii(jj+1);
                % statistics
                result.rp.iterations = result.rp.iterations + iterations;
                
                % handle constraint violations
                result.violations.GAM_DV_violation(jj+1) = GAM_DV_violation;
                result.violations.GAM_Vinf_violation(jj+1) = GAM_Vinf_violation;
                result.violations.rp_violation(jj+1) = rp_violation;                  
                
                % failure
                if failure, result = all_violated(result); return; end  
                
                % set additional DeltaV to zero
                result.additional_DeltaV(jj) = 0;
                
                % set V_arrival, departure
                result.V_arrival(jj+1,:)   = V2;
                result.V_departure(jj+1,:) = output.Vp_vec(2, :); 
                
                % set Vinf_min, plus
                result.Vinf_min(jj,:)  =  output.Vp_vec(1, :); 
                result.Vinf_plus(jj,:) =  output.Vp_vec(2, :);
                
                % DeltaV at central body
                result.DeltaVs(jj) = output.DeltaV;
                % insert location of pericenter
                result.r_arrival(jj,:)   = output.rp_vec;
                result.r_departure(jj+1,:) = output.rp_vec;
            end
            
        elseif exitflag == -1 % problem is impossible
            result.central_body.impossible = 1;
            result = all_violated(result); return;
        else % Solving CBF problem failed
            result.central_body.failure = 1;
            result = all_violated(result); return;
        end
    end % main CBF loop
    
    %% wrapper functions
    
    % wrapper function for general swingby's
    function constraint = CBF_swingby_wrapper(params, V_arrival, V_departure, V_GAMbody, ii)        
        % evaluate GAM function
        [DeltaV, additional_DeltaV, ig,ig,ig,ig, rp_violation,...
            GAM_DV_violation, GAM_Vinf_violation, ig, failure] = swingby(...
            params, V_arrival, V_departure, V_GAMbody, ii);%#ok               
        % value of the constraint function
        constraint = [rp_violation(:); GAM_DV_violation(:); GAM_Vinf_violation(:)];
        if failure, constraint = inf(size(constraint)); end        
    end % CBF swingby wrapper
    
    % wrapper function for the CBF, launch
    function constraint = CBF_launch_wrapper(params, V_departure, VGAM)
        % evaluate launch function  
        % ??? FUTURE WORK: launch cost function
        [ig,ig, constraint] = launch(V_departure, VGAM, ...
            params.max_C3*(1 + params.C3LoverD_tolerance/100));%#ok
    end % CBF launch wrapper
    
    % wrapper function for the CBF, target
    function constraint = CBF_target_wrapper(params, V_arrival, VGAM)%#ok        
        constraint = 0;
        % ??? FUTURE WORK: target cost function
        % Vinf_target = target(V_arrival, VGAM);        
    end % CBF target wrapper 
    
end % Central Body Flyby

% general swingbys
function [DeltaV, additional_DeltaV, rp, Vinfin, Vinfout, turnangle, rp_violation,... 
        GAM_DV_violation, GAM_Vinf_violation, iterations, failure] = swingby(...
        params, V_arrival, V_departure, V_GAMbody, ii)
    
    % initialize
    GAM_DV_violation = 0;       GAM_Vinf_violation = 0;
    failure          = false;   iterations         = 0;
    rp_violation     = 0;       DeltaV             = 0;
    
    % get parameters for the current GAM-body
    mu_GAMbody = params.GMs(ii-1); % GM-value 
    minradius  = params.mean_Radii(ii-1) + params.min_alts(ii-1); % minimum possible [rp]
    VescSOI2   = params.VescSOI(ii-1)^2; % Vesc SOI-error factor
        
    % convert to body-centric velocities
    Vinfout = V_departure - V_GAMbody;
    Vinfin  = V_arrival   - V_GAMbody;

    % initial magnitudes
    Vinfout_dot = Vinfout*Vinfout.';   % dot-products with self
    mVinfout    = sqrt(Vinfout_dot);   % magnitude
    Vinfoutunit = Vinfout/mVinfout;    % unit vector
    Vinfin_dot  = Vinfin*Vinfin.';     % dot-products with self
    mVinfin     = sqrt(Vinfin_dot);    % magnitude
    Vinfinunit  = Vinfin/mVinfin;      % unit vector

    % correct with Vesc at the edge of the SOI
    Vinfout = sqrt(Vinfout_dot + VescSOI2);
    Vinfout = Vinfoutunit * Vinfout;
    Vinfin  = sqrt(Vinfin_dot + VescSOI2);
    Vinfin  = Vinfinunit * Vinfin;

    % re-compute magnitudes
    Vinfout_dot = Vinfout*Vinfout.';  % dot-products with self
    mVinfout    = sqrt(Vinfout_dot);  % magnitudes
    Vinfin_dot  = Vinfin*Vinfin.';    % dot-products with self
    mVinfin     = sqrt(Vinfin_dot);   % magnitudes

    % Calculate turnangle. 
    % Make 100.4% sure it's in (-1 <= dV(+)�dV(-) <= +1)
    turnangle = acos( max(-1,min(1,Vinfin*Vinfout.'/mVinfout/mVinfin)) );
    
    % compute maximum allowable bending angles
    am = mu_GAMbody/Vinfin_dot;    ap = mu_GAMbody/Vinfout_dot;
    maxturnangle = asin(am/(am+minradius)) + asin(ap/(ap+minradius));
    
    % in the general case, there is no additional DeltaV
    additional_DeltaV = 0;

    % if the turnangle is larger than the maximum, and the GAM is
    % executed powered, an additional Delta-V can be used to still
    % accomplish the required turn angle 
    if (turnangle > maxturnangle) && strcmpi(params.types(ii-1), 'powered')
        % additional bending angle required
        Delta_turnangle = turnangle - maxturnangle;                        
        % calculate additional Delta-V required to account for
        % the additional bending
        additional_DeltaV = sqrt(Vinfin_dot + Vinfout_dot ...
            - 2*mVinfin*mVinfout*cos(Delta_turnangle));
        % reset the bad turnangles to the maximum allowed values
        turnangle = maxturnangle;            
        % set the known rp equal to minradius
        rp = minradius;
        
    % otherwise, find the pericenter radius
    else

        % compute the pericenter radius
        % NOTE: using a (1/turnangle)-transformation, the computation takes far
        % less iterations to converge. Moreover, computing (1/turnangle) is
        % usually more stable than computing (turnangle) directly.

        % initial values
        rpplus  = ap/sin(turnangle/2) - ap;
        rpminus = am/sin(turnangle/2) - am;
        rp      = (rpplus + rpminus)/2;
        f       = inf;
        oneoverturnangle = 1/turnangle;

        % do Newton-Raphson iterations
        while (abs(f - oneoverturnangle) > 1e-4)
            % increment iterations
            iterations = iterations + 1;
            % make the procedure more stable
            rp = abs(rp);
            % some substitutions
            apprp = ap + rp;       amprp = am + rp;
            apoamprp = ap/apprp;   amoamprp = am/amprp;
            % output function values
            f  = 1 / (asin(apoamprp) + asin(amoamprp));
            fp = f*f*( ap/apprp^2/sqrt(1 - apoamprp^2) + ...
                am/amprp^2/sqrt(1 - amoamprp^2) );
            % Newton-Raphson step
            rpp = rp; % save old value
            rp  = rp - (f-oneoverturnangle)/fp;
            % try averaging if the number of iterations is getting too high            
            if (iterations == 25), rp = (rpp+rp)/2; end
            % unfortunately, the procedure can still fail. Just set failure
            % to true, and return
            if (iterations > 50), failure = 1; return; end
        end
        
        % [rp] must be larger than the allowed value
        rp_violation = minradius - rp;
        
    end

    % deltaVs at the swingby-body
    DeltaV = sqrt(mu_GAMbody) * abs( sqrt(2/rp+1/am)-sqrt(2/rp+1/ap) );        

    % for unpowered GAMs, the Delta-V MUST be zero. This implies both 
    % Vinf-vectors (in & out) must be equal. If this not the case,
    % reset DeltaV and compute the constraint violation 
    if strcmpi(params.types(ii-1), 'un-powered') && (DeltaV > 0)  
        DeltaV = 0; GAM_Vinf_violation = abs(mVinfin - mVinfout);
    end

    % This deltaV might be larger than the maximum allowable
    % for this particular swingby. In that case, compute the
    % contraint violation
    GAM_DV_violation = DeltaV - ...        
        sqrt(params.max_GAM_DeltaV(ii-1)^2*(1 + params.C3LoverD_tolerance/100));
    
end % normal swingby's

% compute arrival conditions
function [Vinf_target, Vinf_target_violation] = target(V_arrival, V_GAMbody, max_Vinf)
    
    % compute Vinf; initial violation is zero
    Vinf_target = V_arrival - V_GAMbody;  
    
    % compute violation
    Vinf_targetm = sqrt(Vinf_target*Vinf_target.');   % c < 0:
    Vinf_target_violation = Vinf_targetm - max_Vinf;  % negative when constraint is met
    
end % target

% compute launch conditions
function [C3, Vinf_launch, C3_violation] = launch(V_departure, V_GAMbody, max_C3)
    % compute Vinf launch
    Vinf_launch = V_departure - V_GAMbody;
    % compute C3 (= Vinf�Vinf) and violation
    C3 = Vinf_launch*Vinf_launch.';
    C3_violation = C3 - max_C3; % c <=0: negative when constraints are satisfied
end

% set all constraints to maximum violation (for failures etc.)
function result = all_violated(result)
    fields = fieldnames(result.violations);
    for ii = 1:numel(fields)
        result.violations.(fields{ii}) = inf(size(result.violations.(fields{ii}))); 
    end
    result.is_violated = true;      
end % all violated
