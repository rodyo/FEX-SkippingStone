% CALLBACKS()   Contains all callback functions for all controls on the
%               Skipping Stone main window
%
% A Graphical User Interface (GUI), such as Skipping Stone's main window, 
% generally has many so-called User-Interface Controls (UICONTROL()). These
% can be buttons, check boxes, static texts...basically everything that
% makes up a GUI. 
% 
% To tell what a certain UICONTROL should do when a user uses it, the 
% UICONTROL has one or more so-called "callback" functions associated with 
% it. In MATLAB, these can be strings (see EVAL()), built-in MATLAB 
% functions or (of course), manually defined MATLAB functions. 
%
% Although it is common practice to define the callback functions as 
% nested functions or subfunctions in the same M-file, the number of 
% controls on Skipping Stone's main window quickly grew out of control, 
% making the M-file excessively large; it just became too tedious to 
% locate a particular callback function whenever it needed modification. 
%
% That's the purpose of this function -- it is a collection of all 
% callback functions for all UICONTROLS on Skipping Stone's main window. 
% Any particular callback function can conveniently be called by
%
%   [out_arg1, out_arg2, ...] = ...
%       callbacks('callback_function_name', in_arg1, in_arg2, ...)
%
% where 'callback_function_name' is of course the callback-function in 
% question, [in_arg1, "2,...] are the callback-function's input arguments, 
% and [out_arg1, "2,...] the output arguments. 
%
% For a complete list of all callback functions associated with Skipping 
% Stone's main window, type "edit callbacks" at the MATLAB command prompt,
% and scroll down to the beginning of the function's definition. 
%
% See also main, modify_settings.

function varargout = callbacks(funfcn, varargin)
%{
    last edited 19/Dec/2009

    NOTE: This function was written with code-folding of cells, comment 
    blocks and functions in mind. If you have these disabled, it will be 
    quite uncomfortable to read. You can change these options in 

          File > Preferences > Code Folding

    You can fold/unfold individual blocks more quickly with 

          Ctrl + .         (fold)
          Ctrl + Shift + . (unfold)

    NOTE: function names are CaSe SeNsItIvE!

    NOTE: all functions take at least 2 arguments: the calling UICONTROL's 
          handle, and an eventlist. In most of the callback functions below,
          these are hardly ever used directly, so these arguments have been
          included as (varargin) just to make life less confusing.    
%}
  
    %% Initialize
    
    % declare globals
    global MainWin current_tab current_output_tab
    global launch_tab sequence_tab arrival_tab algorithms_tab output_tab
    global Pareto_tab trajectory_tab central_body_speed post_processing 
    global BATCH_optimization optimization_statistics 
        
    % get all appdata
    settings    = getappdata(MainWin, 'settings'   );
    handles     = getappdata(MainWin, 'handles'    );
    environment = getappdata(MainWin, 'environment');
    model       = getappdata(MainWin, 'model'      );
    constants   = getappdata(MainWin, 'constants'  );
    calculation = getappdata(MainWin, 'calculation');
    
    % some abbreviations
    lt = handles.tab(launch_tab);     st  = handles.tab(sequence_tab);
    at = handles.tab(arrival_tab);    ot  = handles.tab(output_tab);     
    ct = handles.tab(current_tab);    alg = handles.tab(algorithms_tab);
    
    % convert to function handle
    % NOTE: STR2FUNC() doesn't work in this context!
    [ig,funfcn] = evalc(['@',funfcn]);%#ok
    
    % CALLBACKS() must have at least ONE output argument in order to use
    % multiple callbacks per UICONTROL(). But since only a few functions
    % actually have output arguments, it's easiest to make exceptions for 
    % just those functions    
    if any(strcmpi(func2str(funfcn), {'callbacks/disable_all';'post_processor'}))
        [varargout{1:nargout}] = feval(funfcn, varargin{:});    
    else
        feval(funfcn, varargin{:});        
        varargout{1} = [];
    end
    
    % write (possibly changed) model,calculation & settings to MainWin's appdata 
    % (NOTE: this construction is pretty convenient, but it doesn't 
    % allow one to use MODIFY_SETTINGS('change single_setting') in any of 
    % the callback functions; settings changed there will get overwritten 
    % by the following statements)
    setappdata(MainWin, 'settings'   , settings);  
    setappdata(MainWin, 'calculation', calculation);
    setappdata(MainWin, 'model'      , model);  
    
    %% Callbacks for menu options
    
    % the about box
    function about(varargin)%#ok
                
        % rename for convenience
        progname     = environment.program_name;
        author       = environment.authors;
        contributors = environment.contributors;
        
        % version info
        Version = num2str(main([]));
        
        % xXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx
        % Create About window below
        % xXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx
        
        % build the author string
        author_string = [];
        for ii = 1:size(author, 1)
            author_string = [author_string;{...
                [' Name         : ', author{ii, 1}];
                [' E-mail       : ', author{ii, 2}];
                [' Affiliation  : ', author{ii, 3} ];
                ['                   ', author{ii, 4} ];
                ['                   ', author{ii, 5} ];
                '                   ';
                '                   '}];%#ok
        end
        
        % build the contributors string
        contributor_string = 'With contributions by: ';
        for ii = 1:length(contributors)
            if ii ~= length(contributors)
                contributor_string = [contributor_string, contributors{ii}, ', '];%#ok
            else
                contributor_string = [contributor_string, 'and ', contributors{ii}, '.'];%#ok
            end
        end
        
        % complete string to be displayed
        about_string = [{['                 ', progname, ' Written by'];
            ' '};
            author_string;
            contributor_string;
            {' ';
            ' ';
            ['Version: ',Version];
            ' ';
            ' ';
            ['           Copyright 2009, all rights reserved. ', progname, ' is '];
	    '            an ongoing project of M.Sc. students at the TU-Delft. ';
 	    '            Currently, it is in the early development phase, and most ';
	    '            likely prone to some bugs. ';
            ' ';
            ' ';
            ' Please report any bugs to Rody.';
            ' ';
            ' ';
            ' ';
            ' ';
            ' ';        % These lines are included to make room for the logo
            ' ';
            ' ';
            ' ';
            ' ';
            ' ';}];
        
        % simply use msgbox
        msgbox(about_string, ['About ', progname]);
        
        % Draw logo        
        % NOTE: the about box will not run when the logo is missing,
        % moved or renamed. Use a TRY-CATCH-END block to prevent this
        try
            % read logo
            logo = imread([environment.pathing.datadir, filesep, 'logo.png']);
            % draw axes
            logo_axes = axes('Position', [.25 .15 .5 .14]);
            % set all whitish-pixels to background color
            bgcolor = environment.colors.window_bgcolor;
            transpinds = (logo(:, :, 1) > 250) & (logo(:, :, 2) > 250) & (logo(:, :, 3) > 250);
            logo1 = logo(:, :, 1); logo1(transpinds) = 255*bgcolor(1);
            logo2 = logo(:, :, 2); logo2(transpinds) = 255*bgcolor(2);
            logo3 = logo(:, :, 3); logo3(transpinds) = 255*bgcolor(3);
            logo = cat(3, logo1, logo2, logo3);
            % draw logo
            imagesc(logo, 'Parent', logo_axes);
            axis off
        catch ME,ME; %#ok<VUNUS>
            warndlg(ME.message, 'Logo not found');
        end % draw logo
        
    end % about
    
    % re-load the default settings
    function load_defaults(varargin)%#ok
        % ask to save the current settings
        if strcmpi(get(MainWin, 'UserData'), 'dirty')
            saveyn = questdlg('Do you want to save the current settings?', ...
                     'Save current settings', 'Yes', 'No', 'Yes');
            if strcmpi(saveyn, 'yes') && save_settings
                set(MainWin, 'UserData', 'clean'), end
        end
        % reload the defaults values
        [environment, model, constants, calculation, settings] = ...
            modify_settings('default_values', ...
            [fileparts(mfilename('fullpath')), filesep, '..']);       
        setappdata(MainWin, 'settings'   , settings   );
        setappdata(MainWin, 'handles'    , handles    );
        setappdata(MainWin, 'environment', environment);
        setappdata(MainWin, 'model'      , model      );
        setappdata(MainWin, 'constants'  , constants  );
        setappdata(MainWin, 'calculation', calculation);
        % and change all settings accordingly
        modify_settings('change_all_settings');
        % hide the panel made by the postprocessor
        if isfield(at.post_processing, 'additional_controls') && ...
           isfield(at.post_processing.additional_controls, 'panel')
            set(at.post_processing.additional_controls.panel(:), 'visible', 'off');            
        end
    end % load defaults
        
    % load settings
    function load_settings(varargin)%#ok
        [settings_file, settings_path] = uigetfile({'*.cfg'}, ...
            ['Load ', environment.program_name, ' configuration settings']);        
        % if cancel was not pressed
        if ischar(settings_file) 
            % check if results are unsaved
            % if so, ask the user to save them
            if strcmpi(get(MainWin, 'UserData'), 'dirty')
                answer = questdlg({'Results from previous optimization have not been saved!';
                    'These will be lost when loading new settings.';
                    'Do you first want to save your results?'}, ...
                    'Unsaved optimization results', 'Yes', 'No', 'Yes');
                if strcmpi(answer, 'yes'), save_results; end
            end              
            % load file into workspace
            try
                new_settings = load([settings_path, settings_file], ...
                    'settings', '-MAT');
            catch ME,ME; %#ok<VUNUS>
                errordlg({'Unable to load file:'; ME.message}, 'Load failed')
                return
            end            
            % check if all parameters were actually saved in this file
            if ~isfield(new_settings, 'settings')
                errordlg({'No program settings found in file.';
                    'Please select a different file.'},...
                    'No settings found');
            else            
                % assign new data
                settings = new_settings.settings;
                % save HERE (otherwise, we can't use change_all_settings)
                setappdata(MainWin, 'settings', settings);
                % re-assign all values to all input fields
                modify_settings('change_all_settings');
            end
        end        
    end % load settings
    
    % save settings
    function success = save_settings(varargin) %#ok<VANUS>
        [settings_file, settings_path] = uiputfile({'*.cfg'}, ...
            ['Save ', environment.program_name, ' configuration settings']);        
        % if cancel was not pressed
        if ischar(settings_file) && ischar(settings_path)   
            % save the settings
            try
                save([settings_path, settings_file], 'settings', '-MAT');
                success = true;
            catch ME,ME; %#ok<VUNUS>
                errordlg({'Unable to save file:'; ME.message}, 'Save failed')
                success = false;
            end         
        else
            success = false;
        end
    end % save settings
    
    % load results
    function load_results(varargin)%#ok
        [results_file, results_path] = uigetfile({'*.res'}, ...
            ['Load ', environment.program_name, ' calculation results']);        
        % if cancel was not pressed
        if ischar(results_file) 
            % check if results are unsaved
            % if so, ask the user to save them
            if strcmpi(get(MainWin, 'UserData'), 'dirty')
                answer = questdlg({'Results from previous optimization have not been saved!';
                    'These will be lost when loading other results.';
                    'Do you first want to save your results?'}, ...
                    'Unsaved optimization results', 'Save', 'Discard', 'Save');
                if strcmpi(answer, 'Save'), save_results; end
            end              
            % load file into workspace
            try
                progress_bar(0, 'Loading results...');
                new_calculation = load([results_path, results_file], ...
                    'calculation', '-MAT');
            catch ME,ME; %#ok<VUNUS>
                errordlg({'Unable to load file:'; ME.message}, 'Load failed')
                return
            end            
            % check if all parameters were actually saved in this file
            if ~isfield(new_calculation, 'calculation')
                progress_bar(0, 'No results found.');
                errordlg({'No optimization results found in file.';
                    'Please select a different file.'},...
                    'No results found');
            else            
                % assign new data
                calculation = new_calculation.calculation;
                setappdata(MainWin, 'calculation', calculation);
                % set flag to clean
                set(MainWin, 'UserData', 'clean');
                % and draw results
                callbacks('showtab', output_tab); % switch to output window
                generate_output('embedded');      % and display results                
            end
            progress_bar('')
        end   
    end % load results
    
    % save results
    function save_results(varargin) %#ok<VANUS>
        % first some checks
        if isempty(calculation.results.first_order.best.solution) && ...
           isempty(calculation.results.second_order.best.solution)
            errordlg({'There are no results to save yet. You should run an ';
                'optimization before saving its results.'}, 'Nothing to save');
            return
        end
        % ask where to save results
        progress_bar(0, 'Saving results...')
        [results_file, results_path] = uiputfile({'*.res'}, ...
            ['Save ', environment.program_name, ' optimization results']);        
        % if cancel was not pressed
        if ischar(results_file) && ischar(results_path)
            % save the settings
            try                
                save([results_path, results_file], 'calculation', '-MAT');                
                set(MainWin, 'UserData', 'clean');                
            catch ME,ME; %#ok<VUNUS>
                errordlg({'Unable to save file:'; ME.message}, 'Save failed')
            end             
        end
        progress_bar('')
    end % save results
    
    % export results
    function export_results(type, varargin)%#ok
        % first some checks
        if isempty(calculation.results.first_order.best.solution) && ...
                isempty(calculation.results.second_order.best.solution)
            errordlg({'There are no results to save yet. You should run an ';
                'optimization before saving its results.'}, 'Nothing to save');
            return
        end
        % ask where to save results
        progress_bar(0, 'Saving results...')
        [results_file, results_path] = uiputfile({['*.' type]}, ...
            ['Export ', environment.program_name, ' optimization results']);
        % if cancel was not pressed
        if ischar(results_file) && ischar(results_path)
            
            % detect type of export
            switch lower(type)
                case {'tab' 'csv'}
                    % set delimiter
                    if strcmp(type, 'tab'), delim = '\t';  %#ok<NASGU>
                    else delim = ',';                      %#ok<NASGU>
                    end
                    % try to DLMwrite "calculation"
                    try
                        % TODO: dlmwrite() can't write structures...
                        % we'll have to call dlmwrite recursively for each
                        % fieldname of [calculation], and prepend the
                        % matrix with a label:
                        error(' '); %#ok<ERTAG>
                        %dlmwrite([results_path,results_file], calculation, ...
                        %    'delimiter', delim);
                    catch ME %#ok<MUCTH,NASGU>
                        % NOTHING YET..
                        not_yet_done;
                    end
                case 'stk'
                    % Sorry, not yet...
                    not_yet_done;
                otherwise
                    % shouldn't ever happen, but oh well...
                    errordlg(['Unknown filetype for export: ''' type '''.'],...
                        'Unknown export type');
            end
        end
    end
    
    % import results
    function import_results(type, varargin)%#ok
        
        [results_file, results_path] = uigetfile({['*.' type]}, ...
            ['Import ',...
            environment.program_name,...
            ' calculation results']);         %#ok<NASGU>
        
        % if cancel was not pressed
        if ischar(results_file) 
            % check if results are unsaved
            % if so, ask the user to save them
            if strcmpi(get(MainWin, 'UserData'), 'dirty')
                answer = questdlg({'Results from previous optimization have not been saved!';
                    'These will be lost when loading other results.';
                    'Do you first want to save your results?'}, ...
                    'Unsaved optimization results', 'Save', 'Discard', 'Save');
                if strcmpi(answer, 'Save'), save_results; end
            end   
            % detect type of import
            switch lower(type)
                case {'tab' ;'csv'}
                    % set delimiter
                    if strcmp(type, 'tab'), delim = '\t';  %#ok<NASGU>
                    else delim = ',';                      %#ok<NASGU>
                    end
                    % try to DLMwrite "calculation"
                    try
                        % TODO: dlmread() can't read structures...
                        % we'll have to call dlmread recursively for each
                        % fieldname of [calculation], depending on the 
                        % prepended label (see export_results):
                        error(' '); %#ok<ERTAG>
                        %dlmread([results_path,results_file], delim);
                    catch ME %#ok<MUCTH,NASGU>
                        % NOTHING YET..
                        not_yet_done;
                    end
                case 'stk'
                    % Sorry, not yet...
                    not_yet_done;
                otherwise
                    % shouldn't ever happen, but oh well...
                    errordlg(['Unknown filetype for import: ''' type '''.'],...
                        'Unknown import type');
            end
        end
    end
    
    %% Callbacks from launch panel
   
   % enable/disable SEP/NEP options
    function change_power_supply(varargin)%#ok
        % enable/disable appropriate controls        
        switch varargin{2}.NewValue
            case lt.power_radio(1)
                % disable NEP-options
                set(lt.NEP, 'enable', 'off');
                % enable SEP-options
                set(lt.SEP, 'enable', 'on');
                % disable or enable the jettison-mass box
                enable_jettison
            case lt.power_radio(2)
                % disable NEP-options
                set(lt.NEP, 'enable', 'on');
                % enable SEP-options
                set(lt.SEP, 'enable', 'off');
                % disable Jettisson
                set(lt.Jettison, 'enable', 'off');
        end
    end
        
    % enable/disable Panel mass edit box
    function enable_jettison(varargin) %#ok<VANUS>
        % enable/disable appropriate controls        
        if get(lt.SEP(7), 'value')
            set(lt.Jettison, 'enable', 'on');
        else
            set(lt.Jettison, 'enable', 'off');
        end      
    end            
       
    % switch engine type
    function engine(varargin)%#ok
        % change the contents and visibility of the edit-boxes        
        switch varargin{2}.NewValue
            case lt.propulsion(1) % High-thrust
                % enable/disaple appropriate engine controls
                set(lt.Sail, 'enable','off'); 
                set(lt.High, 'enable','on'); 
                set(lt.Ion , 'enable','off');
                set(lt.Isp, 'enable','on');                            
                % also change Isp-seting
                set(lt.High(2),...
                    'string'    , settings.propulsion.high_thrust.Isp,...
                    'callback'  , @(varargin) modify_settings('change_single_setting', ...
                        'propulsion.high_thrust.Isp', [], varargin{:}));
                % change the different Isp
                captureIsp;                
            case lt.propulsion(2) % Low-thrust: Ion
                % enable/disaple appropriate controls
                set(lt.Sail,'enable','off');
                set(lt.High,'enable','on');
                set(lt.Ion ,'enable','on');
                set(lt.Isp, 'enable','on');                
                % also change Isp-seting
                set(lt.High(2),...
                    'string'    , settings.propulsion.ion_engine.Isp,...
                    'callback'  , @(varargin) modify_settings('change_single_setting',...
                        'propulsion.ion_engine.Isp', [], varargin{:}));
                % change the different Isp
                captureIsp;                
            case lt.propulsion(3) % Low-thrust: Sail
                set(lt.Sail,'enable','on');
                set(lt.High,'enable','off');
                set(lt.Ion ,'enable','off');
                set(lt.Isp,'enable','off');               
        end             
        % also adjust the controls in the sequence
        set_propulsiondependent_GAM_controls;  
    end
   
    % enable/disable capture Isp
    function captureIsp(varargin) %#ok<VANUS>
        % if checkbox is on        
        if get(lt.Isp(1), 'value')
            set(lt.Isp(2:end), 'enable', 'on')            
        % otherwise
        else
            set(lt.Isp(2:end),'enable', 'off');
            set(lt.Isp(3), 'string', get(lt.High(2), 'string'));
        end
    end
    
    %% Callback from sequence panel
    
    % Check & change dates for the launch window
    function check_dates(varargin)%#ok   
        % rename the value that's changed
        changed_value = varargin{1};    
        % If one of the years OR one of the months has changed:
        % check if month is february. If so, check if year 
        % is leap year. If so, assign 29 days, otherwise 28 days.
        for ii = 1:2
            if any(changed_value == st.launch_window.year(ii) | ...
                   changed_value == st.launch_window.month(ii))
                if get(st.launch_window.month(ii), 'Value') == 2
                    initial_year = get(st.launch_window.year(ii), 'string');
                    initial_year = str2double(initial_year(1, :));
                    year = get(st.launch_window.year(ii), 'Value') + initial_year - 1;
                    if mod(year, 4) == 0
                        set(st.launch_window.day(ii), 'String', 1:29);
                    else
                        set(st.launch_window.day(ii), 'String', 1:28);
                    end
                    return % if february, return
                end
            end
        end   
        % If one of the months has changed:
        % Assign either 30 or 31 days to days-popup, 
        % depending on choice of month        
        for ii = 1:2
            if changed_value == st.launch_window.month(ii)
                % months with 31 days
                if any(get(st.launch_window.month(ii), 'Value') == [1, 3, 5, 7, 8, 10, 12])
                    set(st.launch_window.day(ii), 'String', 1:31);
                    return
                end                
                % months with 30 days
                if any(get(st.launch_window.month(ii), 'Value') == [4, 6, 9, 11])
                    % don't forget to put the cursor on the last position
                    if get(st.launch_window.day(ii) , 'Value') == 31
                        set(st.launch_window.day(ii), 'Value', 30);
                    end
                    % and reset the day
                    set(st.launch_window.day(ii), 'String', 1:30);
                end
            end
        end        
    end  
    
    % enable or disable swingby's according to user's selection
    function switch_bodies(which_GAM_body, varargin) %#ok<VANUS>
        % and what body was selected
        what_body = get(st.GAM.body(which_GAM_body), 'value') - 1;        
        % if not 'none' was selected
        if (what_body ~= 0)            
            % enable all of its controls
            if (which_GAM_body ~= numel(st.GAM.body))
                set(st.GAM.body(which_GAM_body+1), 'enable', 'on');
                % also do recursive call to switch previously 
                % selected bodies back on
                switch_bodies(which_GAM_body+1);
            end
            structfun(@(x) set(x(which_GAM_body), 'enable', 'on'), st.GAM)
            % disable those that are not applicable for the current
            % propulsion or GAM-type settings  
            set_propulsiondependent_GAM_controls;
            change_GAM_type(which_GAM_body);            
            % change the minimum altitude 
            set(st.GAM.minalt(which_GAM_body), 'string', ...
                settings.GAM.min_altitude(what_body, which_GAM_body));
            change_minimum_altitude(which_GAM_body, 'nocheck');
            % change UB/LB to current setting
            set(st.GAM.TOF_LB(which_GAM_body), 'string', ...
                settings.GAM.TOF_LB(what_body+1,which_GAM_body));
            set(st.GAM.TOF_UB(which_GAM_body), 'string', ...
                settings.GAM.TOF_UB(what_body+1,which_GAM_body));            
        % otherwise, disable all following controls
        else  
            structfun(@(x) set(x(which_GAM_body:end), 'enable', 'off'), st.GAM);
            set([st.GAM.body(which_GAM_body), ...
                 st.GAM.number(which_GAM_body)], 'enable', 'on')
        end 
    end
        
    % change the LB/UB on time of flight
    function change_TOF(index, varargin)%#ok
        % get new settings
        new_TOF_LB = get(st.GAM.TOF_LB(index), 'string');
        new_TOF_UB = get(st.GAM.TOF_UB(index), 'string');
        %^ check values
        good = modify_settings('check_value', str2double(new_TOF_LB)) ||...
               modify_settings('check_value', str2double(new_TOF_UB));
       if good
           % currently selected body
           current_body = get(st.GAM.body(index), 'value');
           % and apply the change           
           settings.GAM.TOF_LB(current_body,index) = str2double(new_TOF_LB);
           settings.GAM.TOF_UB(current_body,index) = str2double(new_TOF_UB);           
       end
    end
    
    % change the LB/UB on time of flight
    function change_TOF_target(varargin)%#ok
        % get new settings
        new_TOF_LB = get(st.target.TOF_LB, 'string');
        new_TOF_UB = get(st.target.TOF_UB, 'string');
        % check values
        good = modify_settings('check_value', str2double(new_TOF_LB)) ||...
               modify_settings('check_value', str2double(new_TOF_UB));
       if good
           % currently selected body
           current_body = get(st.target.body, 'value');
           % and apply the change           
           settings.target.TOF_LB(current_body) = str2double(new_TOF_LB);
           settings.target.TOF_UB(current_body) = str2double(new_TOF_UB);           
       end
    end 
    
    % change the minimum altitudes (a neccessary exception to
    % MODIFY_SETTINGS('change_single_setting'))
    function change_minimum_altitude(which_one, check, varargin) %#ok<VANUS>
        % get currently selected body
        current_body = get(st.GAM.body(which_one), 'value') - 1; % remove 'none'
        % get the new value
        new_value = get(st.GAM.minalt(which_one), 'string');
        % check the inserted value (or not)
        if strcmpi(check, 'check')
            good = modify_settings('check_value', str2double(new_value));
        else
            good = true;
        end
        % if it's good, insert it into its proper index
        if good
            settings.GAM.min_altitude(current_body, which_one) = str2double(new_value);   
        % otherwise, reset to the previous value
        else
            set(st.GAM.minalt(which_one), 'string', ...
                settings.GAM.min_altitude(current_body, which_one));
        end
    end
    
    % disable some GAM-controls, depending on 
    % the selected propulsion type
    function set_propulsiondependent_GAM_controls
        % get the index for the last selected body 
        index = sum(cellfun(@(x) strcmpi(x,'on'), ...
            get(st.GAM.body, 'enable')))-1;
        % get the names of the selected bodies
        selected_bodies = get(st.GAM.body(1), 'string');
        values          = get(st.GAM.body, 'value');
        selected_bodies = {selected_bodies{[values{:}].'}}.';%#ok
        % if any of the names of the selected bodies is equal to that of the
        % model's central body, set its GAM-type to "central body flyby"
        centrals = strcmpi(selected_bodies, model.CentralBody{1});
        set(st.GAM.type(centrals), ...
            'string', {'Central body flyby'},...
            'value' , 1);
        % also change the corresponding setting  
        settings.GAM.type(centrals) = {'Central body flyby'};        
        % get the strings for the other bodies
        strings = get(st.GAM.type, 'string');
        % find those that need to be modified
        noncell_inds = ~cellfun(@iscell, strings);        
        strings(noncell_inds) = {strings(noncell_inds)};        
        numel_strings = cellfun(@numel, strings);          
        % if low thrust has been selected    
        if any(~centrals)
            if settings.propulsion.selected(2) || settings.propulsion.selected(3)
                % disable DSM and Delta-V, and max. total DeltaV boxes
                set([st.GAM.DSM, ...
                    st.target.DSM, ...
                    st.GAM.maxDeltaV], 'enable', 'off')
                set([st.MaxTotalDeltaV; 
                     st.GAM.maxDeltaV(:)], ...
                    'string', 0,...
                    'enable', 'off');
                % remove 'powered' from the available types
                % (only those that need to be modified)
                mod_strings = numel_strings == 3;
                set(st.GAM.type(mod_strings & ~centrals), ...
                    'string', {'un-powered';'aerograv'});
                aerogravs = cellfun(@(x)x==3, get(st.GAM.type, 'value'));
                set(st.GAM.type(mod_strings & ~centrals &  aerogravs), 'value', 2);
                set(st.GAM.type(mod_strings & ~centrals & ~aerogravs), 'value', 1); 
            % if high-thrust is selected
            else
                % enable target DSM box
                set(st.target.DSM, 'enable', 'on');
                % enable maximum total Delta-V, and set it to the max.
                set(st.MaxTotalDeltaV, 'enable', 'on');
                current_Isp = settings.propulsion.high_thrust.Isp;
                current_M0  = settings.launch.launch_mass;
                current_Me  = settings.launch.payload_mass;
                max_DV_possible = Tsjiolkovsky(-1, current_Isp, current_M0, current_Me);
                set(st.MaxTotalDeltaV, 'string', max_DV_possible);
                settings.GAM.constraints.max_DV = max_DV_possible;
                % re-insert 'powered' in the available types                
                % (only those that need to be modified)
                mod_strings = numel_strings <= 2;
                set(st.GAM.type(mod_strings & ~centrals), ...
                    'string', {'un-powered';'powered';'aerograv'});                
                aerogravs = cellfun(@(x)x==2, get(st.GAM.type, 'value'));
                set(st.GAM.type(mod_strings & ~centrals &  aerogravs),...
                    'value', 3);
                set(st.GAM.type(mod_strings & ~centrals & ~aerogravs),...
                    'value', 2);
            end
        end   
        % enable DSM and Delta-V boxes
        for ii = 1:index, change_GAM_type(ii); end 
    end
    
    % disable L/D and Max.DV-controls, depending on the 
    % corresponding GAM-type setting
    function change_GAM_type(number, varargin) %#ok<VANUS>
        % get the new setting
        new_type = get(st.GAM.type(number), 'string');
        if ~iscell(new_type), new_type = {new_type}; end
        new_type = new_type{get(st.GAM.type(number), 'value')};
        % change the setting
        settings.GAM.type(number) = {new_type};
        % enable/disable L/D and max. DV edit boxes accordingly
        if strcmpi(new_type, 'un-powered')
            set(st.GAM.LoverD(number), 'enable', 'off');            
            set(st.GAM.maxDeltaV(number), ...
                    'enable', 'off',...
                    'string', 0); % DON'T change the setting
        elseif strcmpi(new_type, 'powered')
            set(st.GAM.DSM(number), 'enable', 'on');
            set(st.GAM.LoverD(number), 'enable', 'off');
            set(st.GAM.maxDeltaV(number), ...
                    'enable', 'on',...
                    'string', settings.GAM.max_DV(number));
        elseif strcmpi(new_type, 'aerograv')
            set(st.GAM.LoverD(number), 'enable', 'on');
            % only for high-thrust
            prop_setting = get(propulsion, 'value');
            if prop_setting{1}
                set(st.GAM.maxDeltaV(number), ...
                    'enable', 'on',...
                    'string', settings.GAM.max_DV(number));
            end
        end
    end
    
    % maximum Delta-V may not exceed what's possible 
    % from Tsjiolkovskii's equation
    function max_total_DV_check(warn_user, varargin)%#ok
        % extract new value
        New_DV = str2double(get(st.MaxTotalDeltaV, 'string'));
        % retreive the current settings for launch- and 
        % payload mass and Isp
        current_Isp = settings.propulsion.high_thrust.Isp;
        current_M0  = settings.launch.launch_mass;
        current_Me  = settings.launch.payload_mass;
        % the value may not be greater than what is allowable by
        % Tsjiolkovskii's equation
        max_DV_possible = Tsjiolkovsky(-1, current_Isp, current_M0, current_Me);
        if (New_DV > max_DV_possible);
            % warn the user (or not)
            if strcmpi(warn_user, 'warn_user')
                warndlg({'This value for the maximum Delta-V is not possible according';
                         'to the selected values for Isp, launch and payload mass.';
                         ' '; ...
                         'The maximum possible value will be inserted instead.'},...
                         'Max. Delta-V not possible')
            end
            % set the maximum value
            set(st.MaxTotalDeltaV, 'string', max_DV_possible);
            % also change the setting
            settings.GAM.constraints.max_DV = max_DV_possible;
        end
    end
    
    % target body might be "Choose minor planet..."
    function target_body(varargin)%#ok
        % get the selected string
        strings = get(st.target.body, 'string');
        % if it is "Choose Minor Planet"
        if strcmpi(strings{get(st.target.body, 'value')},...
                'Choose Minor Planet...')
            % run the function (the same one as the button)
            MP_target_body;
            
            % ??? REMOVE WHEN IMPLEMENTED
            % ??? reset the field (this option is not yet available)
            set(st.target.body, 'value', 1)
            % ??? REMOVE WHEN IMPLEMENTED
            
        else
            % set the proper upper and lower bounds
            settings.target.body = get(st.target.body, 'value');
            set(st.target.TOF_LB, 'string', ...
                settings.target.TOF_LB(settings.target.body));
            set(st.target.TOF_UB, 'string', ...
                settings.target.TOF_UB(settings.target.body));
        end
    end
    
    % choose a minor planet as target    
    function MP_target_body(varargin) %#ok<VANUS>
        
        % nothing yet...
        not_yet_done;
        
        % general function layout:        
        % (1) construct basic GUI        
        % (2) let the user select an MP        
        % (3) adjust the target.body string;
        % (4) adjust the appropriate setting
        % (5) destroy the GUI, but save the settings
        
    end
    
    % switch model
    function model_selection(varargin)%#ok  
        % switch and load the new model
        switch varargin{2}.NewValue
            case st.selected_model(1)  % Solar system
                % load Solar system parameters
                model = Solar_system_parameters(model);                
                % load mission-specific data for the Solar system
                model = user_Solar_system_parameters(model);
                
            case st.selected_model(2)  % Jovian system                 
                % % load Jovian system parameters
                % model = Jovian_system_parameters(model);                
                % % load mission-specific data for the Solar system
                % model = user_Jovian_system_parameters(model);            
                %???not impplemented yet
                not_yet_done;
                set(varargin{1}, 'SelectedObject', varargin{2}.OldValue)                
                
            case st.selected_model(3)  % Julian system
                % % load Julian system parameters
                % model = Julian_system_parameters(model);                
                % % load mission-specific data for the Julian system
                % model = user_Julian_system_parameters(model);
                % % change all settings
                % change_all_settings;
                %???not impplemented yet
                not_yet_done;
                set(varargin{1}, 'SelectedObject', varargin{2}.OldValue)                
                
        end 
        
        % change all settings
        % (NOT YET NEEDED - JUST UNCOMMENT WHEN THE JULIAN/JOVIAN MODELS
        % ARE READY)
        %change_all_settings;
    end
    
    % load the MP-names if the check box is checked.    
    function MPs_check(uncheck, varargin)%#ok
        % if [uncheck] is true, only force the box to be unchecked
        if uncheck
            % disable checkbox
            set(st.selected_model(4), 'value', 0);
            % also adjust setting
            settings.model(4) = false; 
            % and return 
            return
        end
        % otherwise, take action according to current status of the checkbox
        if get(st.selected_model(4), 'Value')
            % enable updatebutton
            set(st.UpdateMPButton, 'enable', 'on')
            % enable MP-target button
            set(st.MPTargetButton, 'enable', 'on')
            % also adjust the targetbody's string  
            string = get(st.target.body, 'string');
            set(st.target.body, 'string', [string; 'Choose Minor Planet...']); 
            % disable all controls while loading
            [objects, states] = disable_all;
            % load the names
            model = minor_planets_parameters(model, environment, constants, false);
            % save the settings
            settings.model(4) = true;
            % reset all controls
            reset_all(objects, states);
        else
            % disable updatebutton
            set(st.UpdateMPButton, 'enable', 'off')
            % disable MP-target button
            set(st.MPTargetButton, 'enable', 'off')
            % reset string
            string = get(st.target.body, 'string');
            set(st.target.body, 'string', string(1:end-1));                
        end
    end
    
    % load the user database upon checking the checkbox 
    function USR_DB_check(varargin)%#ok
        % ??? NOT YET DONE
        not_yet_done;
        set(varargin{1}, 'value', 0);
        return        
        % preliminary function 
%         if get(varargin{1}, 'Value')
%             % enable loadbutton
%             set(st.LoadUserDBButton, 'enable', 'on')
%             % load database
%             % TODO...
%         else
%             % disable loadbutton
%             set(st.LoadUserDBButton, 'enable', 'off')
%         end
    end
      
    % update the MPCORB minor planet database
    function update_minor_planets_data(ask, varargin)%#ok
                                
        % do we need to ask the user for confirmation?
        if ask
            % pop the question
            getastyn = questdlg('Update asteroid data?', 'Update asteroid data',...
                'Yes', 'No', 'Yes');        
            % if no, use the old file and return
            if strcmpi(getastyn, 'No'), return, end
        end
        
        % first disable all inputs
        [objects, states] = disable_all;
        
        % create progress window
        progress_bar(0, 'Updating minor planet model-file, please wait...'); pause(0.25)
        
        % update the MPCORB               
        renamed = false; 
        step = 0;
        try            
            url     = 'http://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT';
            outfile = fullfile(environment.pathing.datadir,...
                               'asteroids', 'MPCORB.DAT');
            
            pause(1), step = 1;            
            progress_bar(0.2, 'Backing up old file...');       
            if (exist(environment.pathing.MP_filename, 'file') == 2)
                success = movefile(environment.pathing.MP_filename, ...
                                   fullfile(environment.pathing.datadir, 'asteroids', 'MPCORB.BAK'),...
                                   'f');
            else
                progress_bar(0.2, 'Old file does not exist.');       
                success = true;
            end
            assert(success, ' '); % <- force TRY block to fail
            
            renamed = true; 
            step = 2;            
            progress_bar(0.5, 'Downloading latest file, please wait...');            
            if verLessThan('MATLAB', '8.4')
                [outfilename, status] = urlwrite(url,outfile);
            else
                status = 1;
                outfilename = websave(outfile, url);
            end      
            
            assert(status==1 && ~isempty(outfilename),...
                   ' '); % <-- force try/catch to fail
            
            progress_bar(1, 'Downloading latest file, please wait...completed.');
            
        catch %#ok
            progress_bar(0, 'FAILED!');
            switch step
                case 0
                    errordlg({'Connection failed!';...
                             'Is an internet connection available?'},...
                             'Conection failed')
                case 1
                    errordlg({'Could not create back-up of old file.'; ....
                             'Check if a file called `MPCORB.DAT` is located in'; ...
                             fullfile(environment.pathing.datadir, 'asteroids')},...
                             'Backup failed')
                case 2
                    errordlg({'Downloading of `MPCORB.DAT` failed!'; ...
                            'Check if MPC updated the location of the data file.'}, ...
                            'Download failed')
            end
            
            % if the last step failed, try to rename the file 
            % back to its old name
            if renamed
                try
                    success = movefile(...
                        [environment.pathing.datadir, filesep, 'asteroids', filesep, 'MPCORB.BAK'],...
                        environment.MP_filename);
                    if ~success
                        error(' '), end %#ok<ERTAG> % <- force TRY block to fail
                    
                catch %#ok
                    environment.pathing.MP_filename = [environment.pathing.datadir, ...
                        filesep, 'asteroids', filesep, 'MPCORB.BAK'];
                    % save immediately                    
                    setappdata(MainWin, 'environment', environment);
                end
            end
            % if we're here, the procedure failed, so return
            progress_bar('');
            % don't forget to enable everything again before returning
            reset_all(objects, states);
            % and return            
            return
            
        end
        
        % update progress bar
        progress_bar(1, 'MP-model sucessfully updated.'); pause(0.5)
        
        % also reload the names
        progress_bar(1, '(Re)loading MP-names...'); pause(0.5) 
        model = minor_planets_parameters(model, environment, constants, true);        
        
        % enable all controls again
        reset_all(objects, states);
        
    end 
    
    % do BATCH-optimization (upon checking the box)  
    function batch_optimization_check(varargin)%#ok
        % enable/disable the button, and perform appropriate actions
        if get(varargin{1}, 'value')
            set(st.BatchOptimizeButton, 'enable','on');
            % disable all sequence panels
            [GAM_objects, GAM_states] = disable_all(st.swingby_panel);
            % save these in the object's UserData
            set(varargin{1}, 'UserData', {GAM_objects, GAM_states});
            % run BATCH-optimization function 
            batch_optimization_wrapper(varargin{:});
        else
            set(st.BatchOptimizeButton, 'enable','off');
            % re-enable all sequence panels
            objstates = get(varargin{1}, 'UserData');
            reset_all(objstates{1}, objstates{2});
        end  
    end % batch_optimization_check
    
    % wrapper function to call the BATCH-optimization window
    function batch_optimization_wrapper(varargin)
        % run the function
        [accept, settings] = batch_optimization;        
        % If "Cancel" was pressed, reset all fields
        if ~accept || ...
                (isfield(settings.BATCH, 'combinations') && isempty(settings.BATCH.combinations))||...
                (~isfield(settings.BATCH, 'combinations'))
            % reset all controls
            set(st.BatchOptimizeButton, 'enable','off');
            set(st.BatchOptimizeCheck,  'value', 0);
            % re-enable all sequence panels
            objstates = get(varargin{1}, 'UserData');
            set(varargin{1}, 'UserData', []);
            reset_all(objstates{1}, objstates{2});
            % don't forget to change the setting
            settings.BATCH.check = false;
        end
    end % batch_optimization_wrapper 
    
    %% Callbacks for arrival & postprocessing tab
    
    %
    function switch_arrival_type(varargin) %#ok
        selected = get(varargin{1}, 'Value');
        set([at.arrival_constraints{:}], 'Visible', 'off');
        set(at.arrival_constraints{selected}, 'Visible', 'on');
        
    end
    
    % enable or disable load / parameters buttons when 
    % user-costfunction is selected
    function arrival_costfunction_check(varargin)%#ok
        if get(varargin{1}, 'value')
            set(at.usr_costfun_arrival(2:end), 'enable', 'on')
        else
            set(at.usr_costfun_arrival(2:end), 'enable', 'off')
        end  
    end
    
    % enable post-processing
    function enable_postprocessing(varargin)%#ok
        if get(varargin{1}, 'value')
            set(at.post_processing.post_processor, 'enable', 'on');
            % enable any additional GUI-components
            if isfield(at.post_processing, 'additional_controls') && ...
               isfield(at.post_processing.additional_controls, 'panel')
                set(at.post_processing.additional_controls.panel(:), 'visible', 'on');
            end            
        else
            set(at.post_processing.post_processor, 'enable', 'off');
            % disable any additional GUI-components
            if isfield(at.post_processing, 'additional_controls') && ...
               isfield(at.post_processing.additional_controls, 'panel')
                set(at.post_processing.additional_controls.panel(:), 'visible', 'off');
            end
        end
    end
    
    % select the post-processor
    function select_postprocessor(varargin)%#ok
        % get the newly selected post processor
        strings = get(varargin{1}, 'string');
        if ~iscell(strings), strings = {strings}; end
        value   = get(varargin{1}, 'value');
        new_selection = strings{value};        
        % take appropriate action
        if strcmpi(new_selection, '(none)')
            % hide any additional GUI-components
            if isfield(at.post_processing, 'additional_controls') && ...
               isfield(at.post_processing.additional_controls, 'panel')
                set(at.post_processing.additional_controls.panel(:), 'visible', 'off');
            end
            % don't forget to change its setting
            settings.postprocessing.post_processor = 1;
            % and we're done
            return
        else
            % find which post-processor's GUI function we have to call
            index = strcmpi({environment.plugin_info.postprocessors(:).name}, ...
                new_selection);
            % save it
            settings.postprocessing.post_processor = index;
            % and call it
            environment.plugin_info.postprocessors(index).GUI_function_handle(value);
        end
    end % select post processor
    
    % Hide/show "find nearby MPs" post-processor threshold controls
    function find_nearby_MPs_threshold(varargin)%#ok
        switch varargin{2}.NewValue
            case at.post_processing.additional_controls.constant_threshold.radio
                structfun(@(x) set(x, 'enable', 'on'), ...
                    at.post_processing.additional_controls.constant_threshold);
                structfun(@(x) set(x, 'enable', 'off'), ...
                    at.post_processing.additional_controls.variable_threshold);
                set(at.post_processing.additional_controls.variable_threshold.radio,...
                    'enable', 'on');
            case at.post_processing.additional_controls.variable_threshold.radio
                structfun(@(x) set(x, 'enable', 'off'), ...
                    at.post_processing.additional_controls.constant_threshold);
                structfun(@(x) set(x, 'enable', 'on'), ...
                    at.post_processing.additional_controls.variable_threshold);
                set(at.post_processing.additional_controls.constant_threshold.radio,...
                    'enable', 'on');
        end
    end % find nearby MPs threshold
    
    %% Callbacks for algorithms tab
    
    % disable certain global optimizers in case of 
    % multi-objective optimization
    function objectives_selection(varargin) %#ok<VANUS>
        % collect statusses
        selected = [get(alg.MaxPayloadObjective  , 'value')
                    get(alg.MinTOFObjective      , 'value')
                    get(alg.OtherObjectives.check, 'value')];
        % max. payload mass should always be selected
        if ~selected(1)
            set(alg.MaxPayloadObjective, 'value', 1);
            settings.optimize.objectives.max_mass = true;
        end                
        % disable all global optimizers but GODLIKE when 
        % multiple objectives are selected       
        if nnz(selected) > 1
            set(alg.glob_opt_group, 'SelectedObject',...
                alg.global_optimizer(1));            
            set(alg.global_optimizer(2:end), 'enable', 'off');          
            settings.optimize.global.optimizer = [true; false; false];
        % enable them if this is not the case
        else
            set(alg.global_optimizer(2:end), 'enable', 'on');
        end
        % enable/disable "Load..." button when "Other objectives" are
        % selected
        if selected(3)
            set(alg.OtherObjectives.button, 'enable', 'on')
        else
            set(alg.OtherObjectives.button, 'enable', 'off')
        end
    end % objectives selection
    
    % set algorithm parameters
    function optimization_parameters(varargin)%#ok
        
        % GODLIKE 
        if settings.optimize.global.optimizer(1)
            % current settings
            s = settings.optimize.global.optimizer_settings{1};
            % pop the dialog  box
            [new_settings, button] = settingsdlg(...
                 s,...
                'title'     , 'GODLIKE parameters',...
                'Separator' , 'GODLIKE',...                  
                'ItersUb'   , s.ItersUb,...
                'ItersLb'   , s.ItersLb,...
                'MaxIters'  , s.MaxIters,...
                'MaxFunEvals', s.MaxFunEvals,...
                'TolCon'    , s.TolCon,...
                {'Display (debug)';'display'}, {'off';'iter'},...
                'NumStreams', {'Auto'; 1;2;3;4;5},...
                'Number of GAs', {1;2;3;4;'none'},...
                'Number of DEs', {1;2;3;4;'none'},...
                'Number of PSOs', {1;2;3;4;'none'},...
                'Number of ASAs', {1;2;3;4;'none'},...
                {'Automatic popsize'; 'AutoPopsize'}, [true true],...
                'popsize', s.popsize,...
                'Separator', ' ',...
                'Separator', 'GA',...
                    'CrossProb', s.GA.CrossProb,...
                    'MutationProb', s.GA.MutationProb,...
                    'Coding', {'Binary';'Decimal'},...
                    'NumBits', s.GA.NumBits,...
                'Separator', 'DE',...
                    'Flb',  s.DE.Flb,...
                    'Fub',  s.DE.Fub,...
                    'CrossConst', s.DE.CrossConst,...
                'Separator', 'ASA',...
                    'T0' , s.ASA.T0,...
                    'ReHeating', s.ASA.ReHeating,...
                'Separator', 'PSO',...
                    'eta1', s.PSO.eta1,...
                    'eta2', s.PSO.eta2,...
                    'eta3', s.PSO.eta3,...
                    'omega', s.PSO.omega,...
                    'NumNeighbors', s.PSO.NumNeighbors,...
                    'NetworkTopology', {'star'; 'ring'; 'fully_connected'});
                % parse settings
                if strcmpi(button, 'OK')
                    % popsize
                    if new_settings.AutoPopsize
                        new_settings.popsize = []; end
                    new_settings = rmfield(new_settings, 'AutoPopsize');
                    % numstreams
                    if strcmpi(new_settings.NumStreams, 'Auto')
                        new_settings.NumStreams = []; end
                    % algorithms
                    new_settings.algorithms = [...
                        repmat({'GA'},  new_settings.NumberOfGAs, 1)
                        repmat({'DE'},  new_settings.NumberOfDEs, 1)
                        repmat({'PSO'}, new_settings.NumberOfPSOs, 1)
                        repmat({'ASA'}, new_settings.NumberOfASAs, 1)];
                     new_settings = rmfield(new_settings, ...
                      {'NumberOfGAs';'NumberOfDEs';'NumberOfPSOs';'NumberOfASAs'});
                    % algorithm specific
                    new_settings.GA.CrossProb    = new_settings.CrossProb;
                    new_settings.GA.MutationProb = new_settings.MutationProb;
                    new_settings.GA.Coding       = new_settings.Coding;
                    new_settings.GA.NumBits      = new_settings.NumBits;
                    new_settings = rmfield(new_settings, ...
                        {'CrossProb';'MutationProb';'Coding';'NumBits'});
                    new_settings.DE.Flb        = new_settings.Flb;
                    new_settings.DE.Fub        = new_settings.Fub;
                    new_settings.DE.CrossConst = new_settings.CrossConst;
                    new_settings = rmfield(new_settings, ...
                        {'Flb';'Fub';'CrossConst'});
                    new_settings.ASA.T0        = new_settings.T0;
                    new_settings.ASA.ReHeating = new_settings.ReHeating;
                    new_settings = rmfield(new_settings, ...
                        {'T0';'ReHeating'});
                    new_settings.PSO.eta1  = new_settings.eta1;
                    new_settings.PSO.eta2  = new_settings.eta2;
                    new_settings.PSO.eta3  = new_settings.eta3;
                    new_settings.PSO.omega = new_settings.omega;
                    new_settings.PSO.NumNeighbors = new_settings.NumNeighbors;
                    new_settings.PSO.NetworkTopology = new_settings.NetworkTopology;
                    new_settings = rmfield(new_settings, ...
                        {'eta1';'eta2';'eta3';'omega';'NumNeighbors';'NetworkTopology'});
                    % assign new options
                    settings.optimize.global.optimizer_settings{1} = new_settings;
                end            
            
        % repeated Nelder Mead
        elseif settings.optimize.global.optimizer(2)
            % current settings
            s = settings.optimize.global.optimizer_settings{2};
            % pop the dialog box
            [new_settings, button] = settingsdlg(...
                 s,...
                'title'      , 'Globalized Nelder-Mead parameters',...
                'Separator'  , 'Repeated Nelder-Mead',... 
                'MaxIters'   , s.MaxIter,...
                'MaxFunEvals', s.MaxFunEvals,...
                'TolCon'     , s.TolCon,...
                {'Display (debug)';'Display'}, {'off';'iter-detailed'},...
                {'Automatic popsize'; 'AutoPopsize'}, [true true],...
                'popsize', s.popsize...
                );
            % parse settings
            if new_settings.AutoPopsize
                new_settings.popsize = []; end
            new_settings = rmfield(new_settings, 'AutoPopsize');
            if strcmpi(button, 'OK')
                settings.optimize.global.optimizer_settings{2} = new_settings; end
                        
        % repeated Quasi-Newton
        elseif settings.optimize.global.optimizer(3)
            % current settings
            s = settings.optimize.global.optimizer_settings{3};
            % pop the dialog box
            [new_settings, button] = settingsdlg(...
                 s,...
                'title'     , 'Globalized Quasi-Newton parameters',...
                'Separator' , 'Repeated Quasi-Newton',... 
                'MaxIters'  , s.MaxIter,...
                'MaxFunEvals', s.MaxFunEvals,...
                'TolCon'    , s.TolCon,...
                {'Display (debug)';'Display'}, {'off';'iter-detailed'},...
                {'Automatic popsize'; 'AutoPopsize'}, [true true],...
                'popsize', s.popsize...
                );
            % parse settings
            if new_settings.AutoPopsize
                new_settings.popsize = []; end
            new_settings = rmfield(new_settings, 'AutoPopsize');
            if strcmpi(button, 'OK')
                settings.optimize.global.optimizer_settings{3} = new_settings; end           
            
        end
    end % optimization parameters
    
    % select costfunctions defined as plugins
    function select_costfcn_plugin(varargin)%#ok
        % create figure
        scz = get(0, 'ScreenSize');         % put the window in the center of the screen
        scx = round(scz(3)/2-200/2);        % (this will usually work fine, except on some
        scy = round(scz(4)/2-400/2);        % multi-monitor setups)
        costfcn_window = figure(...
            'visible' , 'on',...            % hide the GUI while it is being constructed                 
            'position', [scx scy 300 500], ...
            'DockControls', 'off',...       % force it to be non-dockable
            'menubar' ,'none', ...          % menubar is redefined later
            'toolbar' ,'none', ...          % no toolbar (???possible extention)
            'name'    , 'Additional cost functions', ... % window title
            'NumberTitle', 'off',...        % "Figure 123456789:" just looks corny
            'color'   , environment.colors.window_bgcolor,... % use system-default colorscheme
            'defaultuicontrolfontsize', 8); % default font size
        
        % create listbox
        listbox = uicontrol(...
            'style'     , 'list',...
            'background', environment.colors.edit_bgcolor,...
            'units'     , 'normalized',...  
            'parent'    , costfcn_window,...            
            'position'  , [0.01 0.6 0.98 0.39],...
            'string'    , [{'none'}; {environment.plugin_info.costfuns(:).name}],...
            'callback'  , @change_costfcn);
        
        % description
        description_panel = uipanel(...
            'units'     , 'normalized',...  
            'parent'    , costfcn_window,...     
            'position'  , [0.01 0.15 0.98 0.40],...
            'title'     , 'description');            
        description = uicontrol(...
            'style'     , 'text',...            
            'units'     , 'normalized',...  
            'parent'    , description_panel,... 
            'horizontalalignment', 'left',...
            'position'  , [0 0 1 1]);
        
        % "OK" button
        uicontrol(...
            'units'   , 'normalized',...  
            'parent'  , costfcn_window,...  
            'style'   , 'pushbutton',...
            'string'  , 'OK',...
            'position', [0.2 0.05 0.25 0.05],...
            'Callback', @OK);
        
        % "Cancel button"
        uicontrol(...
            'units'   , 'normalized',...  
            'parent'  , costfcn_window,...  
            'style'   , 'pushbutton',...
            'string'  , 'cancel',...
            'position', [0.6 0.05 0.25 0.05],...
            'Callback', {@delete,costfcn_window});
        
        % initialize listbox and description
        index = 1;        
        if ~isempty(settings.optimize.objectives.other.name)            
            names = {environment.plugin_info.costfuns(:).name};
            index = find(strcmpi(settings.optimize.objectives.other.name, names))+1;            
        end
        set(listbox, 'value', index);
        change_costfcn;
        
        % callback for OK-button
        function OK(varargin) %#ok<VANUS>
            % if "none" was selected, de-select the checkbox and its setting
            if (get(listbox, 'value') == 1)
                % reset fields to empty
                settings.optimize.objectives.other.function = [];
                settings.optimize.objectives.other.name = '';
                settings.optimize.objectives.other.use = false;
                settings.optimize.objectives.other.axis_label = '';
                % and de-select the checkbox
                set(alg.OtherObjectives.check, 'value', false);
                objectives_selection;
                
            % otherwise, adjust the appropriate settings 
            else                
                settings.optimize.objectives.other = ...
                    environment.plugin_info.costfuns(get(listbox, 'value')-1);
                settings.optimize.objectives.other.use = true;
            end   
            % set appdata (??? TODO:general construction doesn't work apparently)
            setappdata(MainWin, 'settings', settings);
            % kill window
            delete(costfcn_window);
        end
                  
        % change description when another costfunction is selected
        function change_costfcn(varargin) %#ok<VANUS>
            % find which one was selected
            selected_costfcn = get(listbox, 'value')-1; % subtract "none" entry
            % change the description
            if(selected_costfcn == 0) % "none"
                set(description, 'string', 'Use no additional cost functions.');
            else % "others"
                set(description, 'string',...
                    environment.plugin_info.costfuns(selected_costfcn).description);
            end
        end
        
    end % select costfunction plugin
    
    %% Callbacks for output tab
    
    % Show/hide specific tabs
    function show_output_tab(whichtab, varargin)%#ok
        % hide previous tab        
        set(ot.tab(current_output_tab).button, ...
            'fontweight', 'normal',...
            'value'     , 0); 
        set(ot.tab(current_output_tab).panel, 'visible', 'off');           
        % show the tab and accentuate the button        
        set(ot.tab(whichtab).button, ...
            'fontweight', 'bold', ...
            'value'     , 1);        
        set(ot.tab(whichtab).panel, 'visible', 'on');        
        % set new current tab 
        current_output_tab = whichtab;
        % also set specific axes as current axes
        switch whichtab
            case Pareto_tab
                set(MainWin, 'currentaxes', ot.tab(whichtab).Pareto_pane(1));
            case trajectory_tab
                set(MainWin, 'currentaxes', ot.tab(whichtab).trajectory_pane(1));
            case central_body_speed 
                set(MainWin, 'currentaxes', ot.tab(whichtab).central_body_speed_pane(1));
            case post_processing
                % do nothing                
            case BATCH_optimization
                %set(MainWin, 'currentaxes', ot.tab(whichtab)....);                
            case optimization_statistics
                % do nothing                
        end
    end % show output tab
    
    % reset infopane
    function reset_infopane(which_infopane, varargin)%#ok
        % just set the string of the infopane equal to its UserData
        set(which_infopane, 'string', get(which_infopane, 'UserData'));
    end
    
    % invert the trajectory plot's colors
    function invert_colors(varargin)%#ok
        % rename axes
        trajectory_axes = ot.tab(trajectory_tab).trajectory_pane(1);
        % get all the elements in the plot
        axes_children = get(trajectory_axes, 'children');
        % get all the texts
        text_labels = axes_children(strcmpi(get(axes_children, 'type'),  'text'));
        % set the plot background color to the inverse of what it's now
        current_color = get(trajectory_axes, 'color');
        if all(current_color == [0,0,0])        
                set(trajectory_axes, 'color', 'w')
                set(text_labels, 'color', 'k')
        elseif all(current_color == [1,1,1])
                set(trajectory_axes, 'color', 'k')
                set(text_labels, 'color', 'w')
        end
    end % invert colors
    
    % plot plots in external figure
    function separate_plots(which_one, varargin)%#ok
        switch which_one
            case Pareto_tab
                generate_output('separate_Paretofront')
            case trajectory_tab
                generate_output('separate_trajectory')
            case central_body_speed
                generate_output('separate_central_speed')
        end
    end % separate plots
        
    %% Other callbacks
    
    % Show/hide specific tabs
    function showtab(whichtab, varargin)%#ok  
        
        if current_tab == whichtab
            return; end
        
        % hide previous tab        
        set(handles.tabbutton(current_tab),...
            'fontweight', 'normal',...
            'value'     , 0); 
        set(ct.panel, 'visible', 'off');           
        
        % show the tab and accentuate the button        
        set(handles.tabbutton(whichtab),...
            'fontweight', 'bold', ...
            'value'     , 1);        
        set(handles.tab(whichtab).panel, 'visible', 'on');        
        
        % set new current tab 
        current_tab = whichtab;   
        
    end % show tab
    
    % disable every control
    function [objects, states] = disable_all(varargin)
        % collect all handles
        if (nargin == 0)
            controls = get(ct.panel, 'children');
            objects  = [controls; handles.tabbutton(:); handles.OptimizeButton];
        else
            objects = get(varargin{1}, 'children');
        end
        % initialize         
        states = repmat({'off'}, size(objects, 1), 1);
        % disable everything, while saving the original state
        for ii = 1:numel(states)            
            % UIPANELS and AXES have no 'enable' property
            if any( strcmpi(get(objects(ii), 'type'), {'uicontrol'}))
                % save the state
                if strcmpi(get(objects(ii), 'enable'), 'on')
                    states{ii} = 'on';
                end
                % disable the control
                set(objects(ii), 'enable','off')
                
            % do recursive call for nested UIPANELS
            else
                [nested_objects, nested_states] = disable_all(objects(ii));
                states  = [states; nested_states];  %#ok
                objects = [objects; nested_objects];%#ok
            end
        end
    end % disable all controls
    
    % reset all controls
    function reset_all(objects, states)
        set(objects(strcmpi(states, 'on')), 'enable', 'on');     
    end % reset all controls    
    
end % all callback functions
