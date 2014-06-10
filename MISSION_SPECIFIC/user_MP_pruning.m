function MPs = user_MP_pruning(MPs, constants)
% USR_MP_PRUNING            Mission-specific pre-pruning of the MP-dataset
%
% This function allows you to pre-prune the Minor Planets data-set (MPCORB)
% according to your mission's specifications. When this function is called,
% you will get a variable named [MPs], in which the following fields are
% defined:
%
%
%   MPs.as                      semi-major axes                          [km]
%   MPs.es                      eccentricities                           [-]
%   MPs.is                      inclinations                             [rad]
%   MPs.Omegas                  longitude of ascending node              [rad]
%   MPs.omegas                  argument of perihelion                   [rad]
%   MPs.E0s                     Eccentric anomalies at J2000.0           [rad]
%   MPs.M0s                     Mean anomalies at J2000.0                [rad]
%   MPs.x0s                     Cartesian coordinates at J2000.0         [km and km/s]
%   MPs.theta0s                 true anomalies at J2000.0                [rad]
%   MPs.bs                      semi-minor axes                          [km]
%   MPs.Cs                      location of the elliptic center          [km]
%   MPs.Us                      unit vector towards perihelion           [1]
%   MPs.Vs                      unit vector towards positive [b]         [1]
%   MPs.Ws                      unit normal vector to Us and Vs          [1]
%   MPs.epoch_MJD               epoch in modified Julian dates           [days]
%   MPs.epoch_days_past_J2000   epoch in days past J2000.0               [days]
%   MPs.number_of_observations  number of observations made of the MP    [#]
%   MPs.number_of_oppositions   number of oppositions of the MP          [#]
%   MPs.first_observed          year the MP was first observed           [year]
%   MPs.last_observed           year the MP was last observed            [year]
%   MPs.types                   type of minor planet (see MPCORB format) [-]
%   MPs.names                   name or readible designation number      [-]
%   MPs.numbers                 minor planet number                      [-]
%   MPs.ns                      mean motions                             [rad/s]
%   MPs.Ts                      orbital periods                          [s]
%   MPs.Hs                      angular momentum vector (r x V)          [km2/s]
%   MPs.Hnorms                  norms of angular momentum vectors        [1]
%   MPs.peris                   perihelion radii                         [km]
%   MPs.apos                    aphelion radii                           [km]
%   MPs.number_of_MPs           total number of MPs                      [#]
%
%   plus whatever fields are added / removed in USR_MP_MODEL(). 
%
%
% Note that each of these fields has a number of elements equal to the
% total amount of minor planets in the data set, i.e., 
%
%   numel(MPs.as)    = 399959 = MPs.number_of_MPs,
%   numel(MPs.types) = 399959 = MPs.number_of_MPs, 
%   ...
%
% This function allows you to remove MPs from the default dataset, based 
% on the demands of your mission. It is generally *NOT* a good idea to
% remove the fields themselves. Instead, you should remove *ENTRIES* from
% each field -- removing the field itself will usually cause cryptic
% errors. 
%
% When you want to include the entire dataset unharmed, just leave this
% function empty. 
      
    % Rody :
    % mission specific pruning for 
    % "Trajectory optimization For A Mission to the Solar Bow Shock 
    % and Several Minor Planets".
        
    % I demand that the orbit solution for minor planets that I include in
    % the set must be accurate enough. That means that single-opposition
    % objects or just badly observed minor planets can be ommitted. As a
    % criterion to prune the others out, I demand that the object is
    % observed at least (150 - 10*rp) times, where rp = a*(1-e) = the MP's
    % perihelion. This ensures relatively nearby objects (main asteroid belt,
    % ~3AU) must have around 120 observations, and further removed objects
    % (Kuiper belt, ~30AU) always stay in the set.
    
    % initialize
    persistent pruned 
    if isempty(pruned), pruned = false; end
    
    % have we pruned before?
    if ~pruned
        
        %% THESIS STUFF
        
%         % plot complete & unpruned dataset
%         figure(1), hold on, plot(MPs.peris/150e6, MPs.is*180/pi, 'r.', 'markersize', 1)
%         xlabel('Perihelion [AU]'), ylabel('Inclination [º]')
        
        %% FIRST PRUNING RULE
        
        % find the indices for those MPs which have been 
        % observed (150 - 10*rp) times 
        good_inds = (MPs.number_of_observations >= (150 - 10*MPs.peris/constants.AU));
                
        % prune out those MPs that are observed less times    
        MPs_before_first = MPs.number_of_MPs;% save value
        MPs = rmfield(MPs, 'number_of_MPs'); % remove this one
        field_names = fieldnames(MPs);       % get field names  
        for ii = 1:numel(field_names)        % loop through all fields
            % show progress
            progress_bar(ii/2/numel(field_names), 'Applying first pruning rule...')
            % remove field
            MPs.(field_names{ii}) = MPs.(field_names{ii})(good_inds, :);
        end
         % re-insert the total number
        MPs.number_of_MPs = numel(MPs.as);
        
        % final progress 
        progress_bar(0.5, ['First pruning complete: ',...
            num2str(MPs_before_first - MPs.number_of_MPs), ' MP''s pruned.'])      
        
        %% MORE THESIS STUFF
                
%         figure(1), plot(MPs.peris/150e6, MPs.is*180/pi,...
%             'marker', '.', ...
%             'linestyle', 'none',...
%             'color', [0 0.7 0], ...
%             'markersize', 1)

        %% SECOND PRUNING RULE
        
        % I also assume that the bow shock is considered to be reached when 
        % the spacecraft reaches the 230 AU distance, within a cone with its
        % tip at the Sun, opening angle 40º and pointing to the location of 
        % the stagnation point. In the mission window 2015 - 2040, there are
        % some minor planets that will never come inside this cone, so they can
        % be thrown out

        % show progress
        progress_bar(0.5, 'Applying second pruning rule...')

        % some constants
        peri_Jupiter  = 740573600;           % Jupiter's "current" perihelion [km]
        half_angle    = 45*pi/180;           % half-angle of the cone (45º)
        SBS_direction = [75.4, -7.5]*pi/180; % direction of the bow shock (ecliptic long/lat)    
        th_lb  = SBS_direction(1) - half_angle;
        th_ub  = SBS_direction(1) + half_angle; % limits on the possible coordinates 
        phi_lb = SBS_direction(2) - half_angle; % for the MPs
        phi_ub = SBS_direction(2) + half_angle;

        % this only goes for MP's that can reach outside Jupiter's orbit
        untouched_MPs      = find(MPs.peris <  peri_Jupiter); 
        MPs_second_pruning = find(MPs.peris >= peri_Jupiter); 
                
        % remove this field for now  
        MPs_before_second = MPs.number_of_MPs; % but save value
        MPs = rmfield(MPs, 'number_of_MPs');         

        % get field names  
        field_names = fieldnames(MPs); 

        % set mission lifetime
        t0   = date2days(2015, 1, 1, 0, 0, 0);
        tend = date2days(2040, 1, 1, 0, 0, 0);

        % times at which to get the MP position (1000 positions)    
        num_times = 1000; % ( = about 10 day step)
        times = linspace(t0, tend, num_times);
        
        % prune the MPs 
        % (this is a relatively small loop; only a few hundred MPs)
        good_MP = MPs_second_pruning; % copy; it's easier
        jj = 0;
        for ii = MPs_second_pruning(:).'

            % show progress
            jj = jj + 1;
            if (mod(jj, 25) == 0)
                progress_bar(0.5 + jj/numel(MPs_second_pruning)/2, ...
                    'Applying second pruning rule...');
            end
            
            % extract the orbital elements
            a  = repmat(MPs.as    (ii), numel(times), 1);
            e  = repmat(MPs.es    (ii), numel(times), 1);
            i  = repmat(MPs.is    (ii), numel(times), 1);
            O  = repmat(MPs.Omegas(ii), numel(times), 1);
            o  = repmat(MPs.omegas(ii), numel(times), 1);
            M0 = MPs.M0s(ii);
            n  = MPs.ns (ii);        
            epoch = MPs.epoch_days_past_J2000(ii);

            % first put the MP at my own [t0]
            DeltaT = t0 - epoch;
            M0 = M0 + n*DeltaT*86400; % DON'T FORGET TO CONVERT TO SECONDS!

            % compute all Ms
            Ms = M0 + n*times*86400;  % DON'T FORGET TO CONVERT TO SECONDS!

            % get coordinates
            [coords(:, 1), coords(:, 2), coords(:, 3)] = ...
                kep2cart(a, e, i, O, o, Ms(:), constants.mu_central, 'M');

            % convert to spherical coordinates
            [th, phi] = cart2sph(coords(:, 1), coords(:, 2), coords(:, 3));
            [th, phi] = deal(th, phi);

            % check contraints & prune the MP if it does not 
            % meet the constraints
            constraint1 = phi > phi_lb & phi < phi_ub;
            constraint2 = th  > th_lb  & th  < th_ub ;   
            if ~any(constraint1 & constraint2)
                good_MP = good_MP(good_MP ~= ii);
            end
        end
        
        % include the "untouched" indices 
        MPs_second_pruning = sort([good_MP(:); untouched_MPs(:)]);
        
        % loop through all fields and remove the bad MPs
        for kk = 1:numel(field_names)
            MPs.(field_names{kk}) = MPs.(field_names{kk})(MPs_second_pruning, :);
        end

        % re-insert the total number
        MPs.number_of_MPs = numel(MPs.as);

        % final progress 
        progress_bar(0.5, ['Second pruning complete: ',...
            num2str(MPs_before_second - MPs.number_of_MPs), ' MP''s pruned.'])
                
        % reset progress bar
        progress_bar('');
        
        % NOW we've pruned
        pruned = true;
        
        %% EVEN MORE THESIS STUFF
                
%         % plot pruned dataset on top of the un-pruned version
%         figure(1), plot(MPs.peris/150e6, MPs.is*180/pi, 'b.', 'markersize', 1)
%         legend('Complete set', 'Set after first pruning', 'Set after second pruning')
%         title([num2str(MPs_before_first - MPs.number_of_MPs), ' out of ',...
%                num2str(MPs_before_first),' MP''s pruned.'])
        
    end % if not pruned before 
        
end % user_MP_pruning
