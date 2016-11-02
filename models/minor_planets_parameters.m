function model = minor_planets_parameters(model, environment, constants, reload)
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
% E-mail     : oldenhuis@gmail.com
%

% TODO; why not merge this with minor_planets_model()?

    % Last Edited: 22/Sep/2009
    
    % keep track of whether the MP-names have been loaded
    persistent loaded
    if isempty(loaded)
        loaded = false; end    
        
    % set proper filename
    environment.pathing.MP_filename = fullfile(environment.pathing.datadir, ...
                                               'asteroids',...
                                               'MPCORB.DAT');
    filename = environment.pathing.MP_filename;
    
    % if they have not been loaded before, load them now
    if (~loaded || reload)
        
        % intialize progress bar
        progress_bar(0, 'Reading minor planet names, please wait...'); 
        pause(0.25)
                
        % check if the datafile is present
        if (exist(filename,'file') == 2)
            
            model = minor_planets_model(model,...
                                        constants,...
                                        environment,...
                                        reload);
                    
        % If the file doesn't exist, run update procedure
        else
            % pop the question here
            ButtonPressed = questdlg({'The MPCORB-database was not found.'
                                      'Do you wish to download it now?'}, ...
                                      'MPCORB.DAT not found', 'Yes', 'Cancel', 'Yes');
                                  
            % when affirmative, update the MP-model
            if strcmpi(ButtonPressed, 'yes')
                callbacks('update_minor_planets_data', false);
                
            % otherwise, uncheck the checkbox andadjust the settings
            else
                callbacks('MPs_check', true);
                progress_bar('')
            end
        end % if file exists       
    end % if not loaded 
    
end % minor planets parameters

