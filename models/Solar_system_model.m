function model = Solar_system_model(seq, t0, model, environment, constants)
% SOLAR_SYSTEM_MODEL   Definitions of various parameters for all the planets.


% Authors
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com
% Affiliation: Delft University of Technology
%              Partly performed at ISAS/JAXA

    % Last Edited: 22/Sep/2009
        
    % NOTE - Any sequence always follows the rule:
    % MERCURY-VENUS-EARTH-MARS-JUPITER-SATURN-URANUS-NEPTUNE
    % succeeded by MOON-CERES-PALLAS-VESTA-PLUTO when relevant
       
    % save whether model has been loaded or not
    persistent loaded       
    
    % number of bodies in the Solar system model
    num_bodies = 14;
    
    % there might be user-defined bodies in [seq]; remove them 
    seq = seq(seq <= num_bodies);
        
    % initialize them if this is the first call
    if isempty(loaded), loaded = false(num_bodies, 1); end
    
    % parameters for progress bar    
    init0 = 0.4;
    init  = 0.7;
    increment = 1 - init;
    
    % load only the required textures    
    if any(~loaded(seq)) 
        
        % how many have not been loaded yet?
        which_ones = unique(seq(:));
        which_ones = which_ones(~loaded(which_ones));
        how_many   = nnz(~loaded(which_ones));
        
        % make sure the JPL/DE405 ephemerides are loaded
        % (always needed here to generate the initial x0's)
        model = initialize_ephemerides_generators(seq, ...
                                                  model,...
                                                  environment,...
                                                  'Solar_System');
        
        % initialize waitbar
        progress_bar(0, 'Loading texture maps...');
        
        % load the texture maps
        texturemapfilenames = {
            'sunmap.jpg';
            'mercurymap.jpg';
            'venusmap.jpg';
            {'earthmap1k.jpg', 'earthcloudmap.jpg', 'earthcloudmaptrans.jpg'};
            'marsmap1k.jpg';
            'jupitermap.jpg';
            {'saturnmap.jpg', 'saturnringcolor.jpg', 'saturnringpattern.gif'};
            {'uranusmap.jpg', 'uranusringcolour.jpg','uranusringtrans.gif'};
            'neptunemap.jpg';
            'moonmap1k.jpg';
            [];              % Ceres has no texture
            [];              % Pallas has no texture
            [];              % Vesta has no texture
            'plutomap1k.jpg';
            };
        
        % also load the Sun's texture
        if ~loaded(1)
            which_ones = [1; which_ones];
            how_many = how_many + 1;
        end
        
        % define a default texture for those bodies without texture maps
        default_texture = ones(512,512,3)/2;
        % initialize errormessage
        errmsg = '';
         
        % load selected texture maps
        for ii = 1:how_many            
            % what's the current body?
            which_body = which_ones(ii);
            % update progress bar
            progress_bar(init0+ (init-init0)/how_many*ii, ...
                ['Loading texture maps, ', model.names{which_body}]);
            % some bodies simply don't have a texture; just use a plain
            % gray texture map in that case
            if which_body > size(texturemapfilenames,1) ||...
               isempty(texturemapfilenames{which_body})
                model.textures{which_body, 1} = default_texture;
                continue;
            end            
            % define file path
            filepath = [environment.pathing.datadir,filesep,'planets',filesep,'textures',filesep];            
            % extract correct element
            texture_name = texturemapfilenames{which_body};            
            % Earth, Saturn and Uranus have more than one file
            if iscell(texture_name)
                for jj = 1:length(texture_name)                    
                    % filename
                    filename = [filepath,texture_name{jj}];  
                    try
                        % read image
                        texture = imread(filename);
                        % IMREAD() reverses these arrays for some idiotic reason, so correct
                        if (size(texture, 3) > 1)  % texturemaps are RGB
                            texture = cat(3, flipud(texture(:, :, 1)), ...
                                flipud(texture(:, :, 2)), flipud(texture(:, :, 3)));
                            texture = cat(3, fliplr(texture(:, :, 1)), ...
                                fliplr(texture(:, :, 2)), fliplr(texture(:, :, 3)));
                        else % transparency maps are grayscale
                            texture = rot90(texture,2);
                        end
                    catch ME,ME %#ok<NOPRT>
                        % insert the default texture if this fails
                        texture = default_texture;
                        % append this error message
                        errmsg = char(errmsg, ME.message);
                    end
                    % insert in model structure
                    model.textures{which_body, jj} = texture;                    
                end
            else     
                try 
                % filename
                filename = [filepath, texture_name];                
                % read image                
                texture = imread(filename);                
                % for some idiotic reason, IMREAD() reverses these arrays, so correct
                texture = cat(3, flipud(texture(:, :, 1)), ...
                    flipud(texture(:, :, 2)), flipud(texture(:, :, 3)));   
                texture = cat(3, fliplr(texture(:, :, 1)), ...
                    fliplr(texture(:, :, 2)), fliplr(texture(:, :, 3)));                
                catch ME,ME %#ok<NOPRT>
                    % insert the default texture if this fails
                    texture = default_texture;
                    % append this error message
                    errmsg = char(errmsg, ME.message);
                end
                % insert in params
                model.textures{which_body, 1} = texture;                
            end
        end
        
        % show texture errors (if any)
        if ~isempty(errmsg)
            errordlg(char('Some textures failed to load. The reason(s) were: ',...
                errmsg), 'Some textures failed to load');
        end
        
        % load the rest of the model
        
        % GM-values ( from http://ssd.jpl.nasa.gov )
        progress_bar(init + increment*1/13, 'Loading Solar system model...');
        model.GMs = [...
            132712439940.0       % Sun  
            22032.09             % Mercury 
            324858.63            % Venus  
            398600.44+4902.801   % Earth+Moon
            42828.3              % Mars
            126686511.00+5959.916+3202.739+9887.834+7179.289 % Jupiter+Io+Europa+Ganymede+Callisto
            37931207.80 + 8978.1382 % Saturn+Titan
            5793966.00           % Uranus
            6835107.00+1427.6;   % Neptune+Triton
            4902.801;            % Moon
            63.2;                % Ceres
            14.3;                % Pallas
            17.8;                % Vesta
            6.4e-9*132712439940];% Pluto + Charon (very uncertain)
        
        % initial conditions at given t0
        % (only needed for ephemerides generation by Kepler elements)
        % x0s - Format: [x0 y0 z0 xdot0 ydot0 zdot0]
        progress_bar(init + increment*2/13, 'Loading Solar system model...');                
        % get the initial positions/velocities
        model.x0s = zeros(num_bodies, 6);
        for i = which_ones(:).'
            
            % exclude Sun
            if (i == 1)
                model.x0s(1, :) = zeros(1,6); 
                continue;
            end
            
            % all other bodies
            if environment.have_DE405_MEX
                
                % JPL has different body order; convert
                if i==1, body = 11;
                else body = i-1; end
                
                % JPL expects JD instead of days past J2000
                t0_corrected = t0 + 2451544.5;
                
                model.x0s(i, :) = JPL_DE405_ephemerides(t0_corrected, body);
            else
                model.x0s(i, :) = JPL_DE405_ephemerides(t0, model.states(i, :));
            end
        end
        
        % as, es, is, omegas, Omegas, thetas
        progress_bar(init + increment*3/13, 'Loading Solar system model...');
        el = cart2kep(model.x0s(which_ones, :), constants.mu_central, 'theta');% N-by-6 vector
        % distribute over the different names
        model.as(which_ones,1) = el(:, 1);  model.Omegas(which_ones,1)  = el(:, 4);
        model.es(which_ones,1) = el(:, 2);  model.omegas(which_ones,1)  = el(:, 5);
        model.is(which_ones,1) = el(:, 3);  model.theta0s(which_ones,1) = el(:, 6);
        
        % initial Mean anomalies
        progress_bar(init + increment*4/13, 'Loading Solar system model...');
        model.M0s(which_ones,1) = etheta2M(model.es(which_ones), ...
            model.theta0s(which_ones));
        
        % initial mean motions
        progress_bar(init + increment*5/13, 'Loading Solar system model...');
        goodinds = false(size(model.as));
        goodinds(which_ones) = model.as(which_ones) ~= 0; % prevent division by zero
        model.ns(goodinds) = sqrt(constants.mu_central ./ ...
            model.as(goodinds).^3);
        
        % initial orbital periods
        progress_bar(init + increment*6/13, 'Loading Solar system model...');
        goodinds = false(size(model.as));
        goodinds(which_ones) = model.ns(which_ones) ~= 0; % prevent division by zero
        model.Ts(goodinds,1) = 2*pi ./ model.ns(goodinds);
        
        % initial eccentric anomalies
        progress_bar(init + increment*7/13, 'Loading Solar system model...');
        model.E0s(which_ones,1) = eM2E(model.es(which_ones),...
            model.M0s(which_ones));
        
        % angular momenta
        progress_bar(init + increment*8/13, 'Loading Solar system model...');
        model.Hs(which_ones, :) = cross(model.x0s(which_ones, 1:3), ...
            model.x0s(which_ones, 4:6), 2);
        model.Hnorms(which_ones,1) = sqrt(sum(model.Hs(which_ones, :).^2, 2));
        
        % apses        
        progress_bar(init + increment*9/13, 'Loading Solar system model...');
        model.peris(which_ones,1) = model.as(which_ones).*...
            (1 - model.es(which_ones));
        model.apos(which_ones,1)  = model.as(which_ones).*...
            (1 + model.es(which_ones));
        
        % Albedo's (Bond albedo's, from wikipedia.org )
        progress_bar(init + increment*10/13, 'Loading Solar system model...');
        model.epsilons = [0; 0.119; 0.65; 0.367; 0.150; 0.343; 0.342; ...
            0.300; 0.290; 0.120];
        
        % spheres of influence
        progress_bar(init + increment*11/13, 'Loading Solar system model...');
        model.SOIs(which_ones,1) = model.as(which_ones) .* ...
            (model.GMs(which_ones)/constants.mu_central).^(2/5);
        
        % Vesc at the edge of the SOI
        progress_bar(init + increment*12/13, 'Loading Solar system model...');
        model.VescSOI(which_ones,1) = ...
            sqrt(2*model.GMs(which_ones) ./ model.SOIs(which_ones));
        
        % all done
        progress_bar(1, 'Solar system model successfully loaded.');
        pause(0.5), progress_bar('');
        
        % model for these bodies is now loaded into memory
        loaded(which_ones) = true;
    end
    
end
