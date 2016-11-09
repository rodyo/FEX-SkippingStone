function model = initialize_ephemerides_generators(seq, model, environment, type)
% initialize all the ephemerides generators   
    
    % XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx
    % INITIALIZE
    % XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx
    
    % save whether model has been loaded or not
    persistent loaded loaded_type
            
    % take care to make seq a col-vector
    seq = seq(:);      
    
    % determine which data to load
    if isempty(loaded) || any(~loaded(seq(seq < numel(loaded))))
        if strcmpi(type, 'Solar_System')
            
            % number of bodies in the Solar system model
            num_bodies = 14;
            
            % stored high-accuracy ephemerides for the planets
            % NOTE: although superfluous, it's a ot easier later on to
            % also load the Sun's statevectors (which are all zero of course).
            statevecData = {[]; @statemercury; @statevenus; 
                @stateearth; @statemars; @statejupiter; @statesaturn;
                @stateuranus; @stateneptune; @statemoon; @stateceres;
                @statepallas; @statevesta; @statepluto};
            
            % store the type of model loaded
            loaded_type = 'Solar_System';
            % remove the Sun (it's position is known)
            seq(seq == 1) = [];
            
        elseif strcmpi(type, 'Jovian_System')
            %??? to be implemented
            
                        
        elseif strcmpi(type, 'Julian_System')
            %??? to be implemented
            
            
            
        end
        
        % see if any user-defined bodies are included in the sequence
        if any(seq > num_bodies)
            
            % find which ones they are
            user_bodies = seq(seq > num_bodies);
            
            % if they are static bodies, we don't need to load ephemerides            
            static_user_bodies = model.static{1}(user_bodies);            
            
            % if there are any non-static bodies, look for an ephemerides file
            % and insert it in the array of ephemerides functions
            if any(~static_user_bodies) 
                % create default function names
                % (default name is "state(body name).m", where the [body
                % name] is extracted (and spaces removed) from [model.names])
                user_names = model.names(user_bodies(~static_user_bodies));
                eph_files  = cellfun(@(x) x(~isspace(x)), user_names, 'uniformoutput', false);
                eph_files  = cellfun(@(x) ['state', x], eph_files, 'uniformoutput', false);
                eph_files_exist = cellfun(@exist, eph_files);
                
                % insert the ones that have been found
                good_inds = user_bodies(~static_user_bodies & eph_files_exist);
                if any(good_inds)
                    [dummy, statevecData{good_inds}] = evalc(['@',eph_files(good_inds)]); %#ok<ASGLU>
                end
                
                % if it doesn't exist, try to generate it using Kepler
                if any(~good_inds)
                    %??? TODO
                    % if that also fails, give an error
                end
                
                % update the number of bodies
                num_bodies = max(user_bodies);
            end    
            
            % if there are any static bodies, just remove their indices from [seq] 
            % (the ephemerides are just returned from [model.static{2}])
            if any(static_user_bodies) 
                seq(seq == user_bodies(static_user_bodies)) = [];                
            end            
        end % if 
        
        % initialize persistent variable
        loaded = false(num_bodies, 1);
        
    end % if (isempty(loaded))
    
    
    % XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx
    % LOAD THE MODEL 
    % XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx
       
    if ~environment.have_DE405_MEX
    
        % load only the required statevecs
        if any(~loaded(seq(seq < numel(loaded))))             
            
            % how many have not been loaded yet?
            which_ones = unique(seq);
            how_many   = nnz(~loaded(which_ones));
            
            % initialize progress bar
            progress_bar(0, 'Interpolating state-vectors, ');
            
            %step = 1/how_many;
            % load selected statevectors
            for ii = 1:how_many
                % what's the current body?
                which_body = which_ones(ii);
                % update progress bar
    %             progress_bar(ii*step, ...
    %                 ['Interpolating state-vectors, ', model.names{which_body}]);
                % load correct datafile            
                stateBody = feval(statevecData{which_body});
                % initialize the number of points
                % (different per datafile, beware)
                numpts = 10*length(stateBody) - 1;
                pts    = 0:10:numpts; % always 10-day step
                % interpolate all 6 elements of the statevector
                for jj = 1:6
                    % interpolate
                    model.states{which_body, jj} = spline(pts, stateBody(:, jj));
                    % only store the breakpoints and coefficients
                    model.states{which_body, jj} = rmfield(model.states{which_body, jj},...
                        {'form', 'pieces', 'order', 'dim'});
                end            
            end
            
            % statevect are now loaded
            loaded(which_ones) = true;

            return;

        % subsequent calls
        else
            % check if the model changed
            if ~strcmpi(loaded_type, type);
                
                % clear the persistent variables
                loaded = []; loaded_type = [];
                
                % and do a recursive call
                model = initialize_ephemerides_generators(seq, model, environment, type);        
                
            end
        end  

    else
        model.states = [];
    end
    
end
