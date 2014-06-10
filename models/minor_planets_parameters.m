function model = minor_planets_parameters(model, environment, reload)
%MINOR_PLANET_PARAMETERS    Definition file for various parameters for asteroids.
%
%   MINOR_PLANET_PARAMETERS only loads the names of
%   all minor planets, since only the names are
%   needed while setting the parameters for the final
%   optimization. Only when the OPTIMIZE-button is
%   pressed will the full model be loaded. 

% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology
%

    % Last Edited: 22/Sep/2009
    
    % keep track of whether the MP-names have been loaded
    persistent loaded
    if isempty(loaded), loaded = false; end    
    
    % set proper filename
    environment.pathing.MP_filename = ...
        [environment.pathing.datadir, filesep, 'asteroids', filesep, 'MPCORB.DAT'];
    
    % if they have not been loaded before, load them now
    if (~loaded || reload)
        
        % intialize progress bar
        progress_bar(0, 'Reading minor planet names, please wait...'); pause(0.25)
                
        % check if the datafile is present
        if (exist(environment.pathing.MP_filename, 'file') == 2)
        
            % open the file
            fid = fopen(environment.pathing.MP_filename); 

            % read blocks of 10.000 at a time to prevent out-of-memory errors
            MP_names = []; counter = 0; num_MPs = 1e4;
            total = 400000; % <- first estimate on the total amount of Minor Planets
            while ~feof(fid)
                if (counter == 0)
                    % skip the header in the first read
                    names = textscan(fid, '%*166s%27s%*9s', num_MPs, ...
                        'headerlines', 39, 'Whitespace', '');
                else
                    names = textscan(fid, '%*166s%27s%*9s', num_MPs, 'Whitespace', '');
                end
                % increase counter
                counter = counter + num_MPs;
                % remove spaces 
                names = strtrim(names);
                % show progress            
                progress_bar(counter/total, [num2str(counter), ' sucessfully read...']);
                % append names to output
                MP_names = [MP_names; names{:}]; %#ok
            end

            % insert names into model
            model.MPs.MP_names = MP_names;

            % flose the file
            fclose(fid);

            % set number of asteroids
            model.MPs.number_of_MPs = numel(model.MPs.MP_names);

            % update and close waitbar
            progress_bar(1, [num2str(model.MPs.number_of_MPs), ' names loaded into memory.']);
            pause(1) % total number of asteroids should be visible for at least SOME time

            % keep them in the workspace until the program is closed,
            % or the datafile is updated
            loaded = true;

            % reset the progress bar
            progress_bar('');
        
        % if the file doesn't exist, run update procedure
        else
            % pop the question here
            ButtonPressed = questdlg({'The MPCORB-database was not found.'
                'Do you wish to download it now?'}, 'MPCORB.DAT not found', 'Yes', 'Cancel', 'Yes');
            % when affirmative, update the MP-model
            if strcmpi(ButtonPressed, 'yes')
                callbacks('update_minor_planets_data', false);
            % otherwise, uncheck the checkbox andadjust the settings
            else
                callbacks('MPs_check', true);
                % also reset progress bar
                progress_bar('')
            end
        end % if file exists       
    end % if not loaded      
end % minor planets parameters
