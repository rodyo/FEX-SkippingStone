% method by Sims and Flanagan
function SimsFlanagan(rs, tfs, Vmins, Vpluss, options, MGA_data, SC_data, muC)

    %% INITIALIZE

    % parse input


    % initial estimate
    %
    % FORMAT:
    %   [t0;
    %    V1(1); DeltaV1(1); DeltaV1(2); ... ; DeltaV1(N1); V2(1); tf(1);
    %    V1(2); DeltaV2(1); DeltaV2(2); ... ; DeltaV2(N2); V2(2); tf(2);
    %    ...
    %    V1(M); DeltaVM(1); DeltaVM(2); ... ; DeltaVM(NM); V2(M); tf(M)];
    %
    % Where



    % objective/constraint functions
%     objective_function  = @(DeltaVs) LT_leg_objective(DeltaVs, SC_data);
%     constraint_function = @(DeltaVs) LT_leg_constraints(r1, r2, tf, DeltaVs, SC_data, muC);

    %% ROUTINE

    % loop through all [rs]
    for i = 1:size(rs,1)-1

        % initialize
        r1 = rs(i,:);
        r2 = rs(i+1,:);

N = ceil( diff(tfs)/options.timestep );
N(mod(N,2)~=0) = N(mod(N,2)~=0)+1;
SC_data.Isp = 3000;
SC_data.P0 = 10e3;
SC_data.m_start = 3000;
SC_data.m_end = 2500;
tf = tfs(2);


        % initial estimates
        Vs = [linspace( Vpluss(i, 1), Vmins(i+1, 1), N+1).',...
              linspace( Vpluss(i, 2), Vmins(i+1, 2), N+1).',...
              linspace( Vpluss(i, 3), Vmins(i+1, 3), N+1).'];
        initial_estimate = diff(Vs)/5;

        % include terminal velocities
        initial_estimate(1,:) = Vs(1, :);
        initial_estimate(end,:) = Vs(end, :);

objective_function  = @(DeltaVs) LT_leg([], [], [], DeltaVs, SC_data, muC, 'objective');
constraint_function = @(DeltaVs) LT_leg(r1, r2, tf, DeltaVs, SC_data, muC, 'constraint');


        % FMINCON options
        options = optimset(...
            'display', 'iter-detailed',...
            'TolCon', 1e-6,...
            'TolX'  , 1e-4,...
            'TolFun', 1e-6,...
            'MaxFunEvals', 1e5,...
            'MaxIter', 1500,...
            'algorithm', 'interior-point');

       % FMINCON call
       [sol, fval, exitflag, output] = fmincon(...
           objective_function, initial_estimate, [],[],[],[],[],[], constraint_function, options);

    end

    % objective function / nonlinear constraint function
    function varargout = LT_leg(r1, r2, tf, DeltaVs, SC_data, muC, mode)

        %% initialize

        % split DeltaVs into terminal velocities and forward and backward part
        V1 = DeltaVs(1, :);     DeltaVs_fwd = DeltaVs(2:N/2, :);
        V2 = DeltaVs(end, :);   DeltaVs_bwd = DeltaVs((end-1):-1:(N/2+1), :);

        % compute end mass
        M_end = tsjiolkovsky(sum(sqrt(sum(DeltaVs(2:end-1,:).^2,2))), ...
            SC_data.Isp, SC_data.m_start, []);

        % cost is simply negative mass
        if strcmpi(mode, 'objective'), varargout{1} = -M_end; return, end

        % constants & initial values
        c  = zeros(sum(N),1); % inequality constraints
        g0 = 9.80665;         % std. grav. accel. at sealevel
        AU = 1.49597870691e8; % 1 AU


        %% DEBUG

if (rand < 0.01), plotit = true;
else plotit = false; end
if plotit, figure(1), clf, hold on, end


        %% progress start and endpoints to matchpoints

        % initialize backwards progression
        ts = linspace(0, tf, N).';              step = ts(2)-ts(1);
        states_fwd = [r1, V1];                  states_bwd = [r2, V2];
        M_fwd = 3000;                           M_bwd = M_end;

        % progress both to matchpoint
        for ii = 1:(N/2-1)

            %% compute inequality constraints

            % calculate distances to Sun (in AU)
            r_fwd = sqrt(states_fwd(1:3) * states_fwd(1:3).')/AU;
            r_bwd = sqrt(states_bwd(1:3) * states_bwd(1:3).')/AU;
            % compute magnitudes of DeltaV
            DeltaVm_fwd = sqrt(sum(DeltaVs_fwd(ii,:).^2,2));
            DeltaVm_bwd = sqrt(sum(DeltaVs_bwd(ii,:).^2,2));
            % available power
            P_fwd = SC_data.P0 / r_fwd/r_fwd;
            P_bwd = SC_data.P0 / r_bwd/r_bwd;
            % max. acceleration possible
            DeltaV_max_fwd = 2*P_fwd*step/M_fwd/g0/SC_data.Isp * 86.4; % (= 60*60*24/1000,
            DeltaV_max_bwd = 2*P_bwd*step/M_bwd/g0/SC_data.Isp * 86.4; %  for unit conversion)
            % inequality constraints
            c(ii    ) = DeltaVm_fwd - DeltaV_max_fwd;
            c(end-ii) = DeltaVm_bwd - DeltaV_max_bwd;

            %% compute conditions at new points

            % new statevector
            states_fwd = progress_orbit(+step, states_fwd, muC);
            states_bwd = progress_orbit(-step, states_bwd, muC);
            % Add DeltaV
            states_fwd(4:6) = states_fwd(4:6) + DeltaVs_fwd(ii, :);
            states_bwd(4:6) = states_bwd(4:6) - DeltaVs_bwd(ii, :);
            % calculate new masses
            M_fwd = tsjiolkovsky(DeltaVm_fwd, SC_data.Isp, M_fwd, []);
            M_bwd = tsjiolkovsky(DeltaVm_bwd, SC_data.Isp, [], M_bwd);


if plotit
plot(states_fwd(:, 1)/150e6, states_fwd(:, 2)/150e6, 'r.')
plot(states_bwd(:, 1)/150e6, states_bwd(:, 2)/150e6, 'b.')
line([states_bwd(:, 1)/150e6, states_bwd(:, 1)/150e6+DeltaVs_bwd(ii, 1)],...
[states_bwd(:, 2)/150e6, states_bwd(:, 2)/150e6+DeltaVs_bwd(ii, 2)])
line([states_fwd(:, 1)/150e6, states_fwd(:, 1)/150e6+DeltaVs_fwd(ii, 1)],...
[states_fwd(:, 2)/150e6, states_fwd(:, 2)/150e6+DeltaVs_fwd(ii, 2)],...
'color','r')
end

        end

        % rename for clarity
        states_matchpoint_fwd = states_fwd;
        states_matchpoint_bwd = states_bwd;

if plotit
    axis equal
    drawnow
end

        %% nonlinear equality constraints

        varargout{1} = c;
        varargout{2} = [...
            M_fwd - M_bwd                % masses at matchpoint should be equal
            DeltaVs_fwd(end,:).'- ...
            DeltaVs_bwd(end,:).'              % Delta-V's at matchpoint must be equal
            states_matchpoint_bwd(:)/1e5 - ...
            states_matchpoint_fwd(:)/1e5]; % positions & velocities should be equal

        %% nonlinear inequality constraints

    end


end
