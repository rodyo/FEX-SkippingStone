% BATCH-optimization GUI
%
%
%
%
function [accept, settings] = batch_optimization(varargin)

    %% Initialize
    
    % get relevant globals
    global MainWin
    
    % get relevant data
    environment = getappdata(MainWin, 'environment');
    model       = getappdata(MainWin, 'model'      );
    settings    = getappdata(MainWin, 'settings'   );
    
    % initialize BATCH-data
    if isfield(settings.BATCH, 'best_of')
        BATCH = settings.BATCH;
    else
        BATCH.combinations = {};
        BATCH.exclusions   = false;
        BATCH.num_combinations = 0;
        BATCH.best_of      = 1;
        BATCH.min_GAMs     = 1;
        BATCH.max_GAMs     = 1;
        BATCH.bodies       = {};
        BATCH.check        = true;
    end
    
    % initially, don't accept settings
    accept = false;    
    % selectable bodies
    body_names = model.names(model.GAMable);    
    % Abbreviated names
    abbrv_names = abbreviate_names;
    
    %% Construct gui
    
    % get screensize
    scz = get(0, 'ScreenSize');                    % put the window in the center of the screen
    scx = round(scz(3)/2-500/2);                   % (this will usually work fine, except on some  
    scy = round(scz(4)/2-600/2);                   % multi-monitor setups) 

    % construct basic GUI
    batch_window = figure;
    set(batch_window,...
        'position', [scx scy 500 600],...          % should be fine
        'renderer', 'zbuffer', ...                 % safest & fastest choice 
        'name'    , 'BATCH-optimization',...       % "Figure 2" is just corny
        'numbertitle', 'off',...
        'visible' , 'off',...                      % don't show until the window is fully drawn
        'units'   , 'normalized',...               % better for resizing etc.
        'resize'  , 'on', ...                      % but just keep it un-resizable for now
        'menubar' , 'none', ...                    % menubar is redefined later
        'toolbar' , 'none', ...                    % no toolbar (???possible extention)
        'windowstyle' , 'modal',...
        'DockControls', 'off',...                  % force it to be non-dockable
        'CloseRequestFcn', @cancel_button,...      % closing the window equals pressing "Cancel"
        'color'   , environment.colors.window_bgcolor,...  % use system-default colorscheme            
        'defaultuicontrolfontsize', 8);            % force font size
        
    % static text    
    help_string = {['For many missions you will quickly find yourself optimizing tens or hundreds ',...
        'of different sequences manually. This process can be automized by selecting all ',...
        ' swingby-bodies to try in the window below. When optimizing with this option enabled, ',...
        environment.program_name, ' will automatically test all combinations possible with those ',...
        'bodies.']};
    explanatory_text = uicontrol(...
        'units'   , 'normalized',...        
        'style'   , 'text',...
        'horizontalalignment', 'left',...
        'position', [0.02 0.89 0.96 0.10]);    
    set(explanatory_text, 'string', textwrap(explanatory_text, help_string));
    
    % list of selectable GAM-bodies
    available_bodies = uipanel(...        
        'title'   , 'Available GAM-bodies',...
        'position', [0.02, 0.32, 0.40, 0.56]);
    available_bodies_list = uicontrol(...     
        'parent'    , available_bodies,...
        'style'     , 'listbox',...
        'units'     , 'normalized',...
        'Background', environment.colors.edit_bgcolor,...
        'string'    , body_names,...
        'position'  , [0.02, 0.02, 0.96, 0.96]);
    
    % list of selected GAM-bodies
    selected_bodies = uipanel(...        
        'title'   , 'Selected GAM-bodies',...
        'position', [0.58, 0.32, 0.40, 0.56]);
    selected_bodies_list = uicontrol(...     
        'parent'    , selected_bodies,...
        'style'     , 'listbox',...
        'units'     , 'normalized',...
        'Background', environment.colors.edit_bgcolor,...
        'string'    , {},...
        'position'  , [0.02, 0.02, 0.96, 0.96]);
    
    % Show number of combinations
    combinations_text(1) = uicontrol(...
        'style'     , 'text',...
        'position'  , [0.20 0.28 0.55 0.02],...
        'units'     , 'normalized',...
        'horizontalalignment', 'left',...
        'string'    , 'Total number of optimizations to perform:',...        
        'enable'    , 'off',...
        'fontweight', 'bold');
    combinations_text(2) =  uicontrol(...
        'style'     , 'text',...
        'position'  , [0.68 0.28 0.5 0.02],...
        'units'     , 'normalized',...
        'horizontalalignment', 'left',...
        'enable'    , 'off',...
        'string'    , '0',...
        'fontweight', 'bold',...
        'ForegroundColor' , 'r');  

    % "Add" button
    uicontrol(...
        'style'   , 'pushbutton',...
        'string'  , 'add >>',...
        'units'   , 'normalized',...
        'position', [0.44 0.65 0.115 0.04],...
        'callback', @add_body);    
    % "Remove" button
    remove_button = uicontrol(...
        'style'   , 'pushbutton',...
        'string'  , '<< remove',...
        'enable'  , 'off',...
        'units'   , 'normalized',...
        'position', [0.44 0.60 0.115 0.04],...
        'callback', @remove_body);
    % "Clear" button
    clear_button = uicontrol(...
        'style'   , 'pushbutton',...
        'string'  , 'Clear',...
        'enable'  , 'off',...
        'units'   , 'normalized',...
        'position', [0.44 0.40 0.115 0.04], ...
        'callback', @clear_all_bodies);

    % Options panel
    options_panel = uipanel(...
        'position', [0.02 0.1 0.96, 0.15]);
        % Use best of [X] optimizations per combination
        uicontrol(...
            'parent'  , options_panel,...
            'units'   , 'normalized',...
            'style'   , 'text',...
            'string'  , 'Use best of',...
            'horizontalalignment', 'left',...            
            'position', [0.02 0.65 0.15, 0.25]);
        best_of = uicontrol(...
            'parent'  , options_panel,...
            'units'   , 'normalized',...
            'style'   , 'edit',...
            'string'  , 1,...
            'Background', environment.colors.edit_bgcolor,...            
            'position', [0.18 0.7 0.15, 0.25],...
            'callback', @show_number_of_optimizations);
        uicontrol(...
            'parent'  , options_panel,...
            'units'   , 'normalized',...
            'style'   , 'text',...
            'string'  , 'optimizations per combination',...
            'horizontalalignment', 'left',...            
            'position', [0.35 0.65 0.45, 0.25]);
        % set maximum/minimum amount of swingby's
        uicontrol(...
            'parent'   , options_panel,...
            'units'    , 'normalized',...
            'style'    , 'text',...
            'position' , [0.02 .23 .23 .25],...
            'horizontalalignment', 'left',...
            'string'   , 'Use minimum of');
        min_GAMs_popup = uicontrol(...
            'parent'  , options_panel,...
            'units'   , 'normalized',...
            'style'   , 'popup',...
            'position', [0.2 .25 .15 .25],...
            'string'  , 0,...
            'Background', environment.colors.edit_bgcolor,...
            'enable'  , 'off',...
            'callback', @show_number_of_optimizations);
        uicontrol(...
            'parent'   , options_panel,...
            'units'    , 'normalized',...
            'style'    , 'text',...
            'position' , [0.36 .23 .20 .25],...
            'horizontalalignment', 'left',...
            'string'   , 'and maximum of');
        max_GAMs_popup = uicontrol(...
            'parent'  , options_panel,...
            'units'   , 'normalized',...
            'style'   , 'popup',...
            'position', [0.55 .25 .15 .25],...
            'string'  , 0,...
            'Background', environment.colors.edit_bgcolor,...
            'enable'  , 'off',...
            'callback', @show_number_of_optimizations);
        uicontrol(...
            'parent'   , options_panel,...
            'units'    , 'normalized',...
            'style'    , 'text',...
            'position' , [0.72 .23 .20 .25],...
            'horizontalalignment', 'left',...
            'string'   , 'swingby''s');
        
    % "exclude combinations" button  
    exclusion_list_button = uicontrol(...
        'style'   , 'togglebutton',...
        'string'  , 'Exclude combinations >>',...
        'enable'  , 'off',...
        'units'   , 'normalized',...
        'position', [0.01 0.02 0.30, 0.05],...
        'callback', @(varargin) exclusion_list(BATCH, varargin{:})); 
    
    % "Cancel" button
    uicontrol(...
        'style'   , 'pushbutton',...
        'string'  , 'Cancel',...
        'units'   , 'normalized',...
        'position', [0.62 0.02 0.15, 0.05],...
        'callback', @cancel_button);
    
    % "OK" button
    OK_button = uicontrol(...
        'style'   , 'pushbutton',...
        'string'  , 'OK',...
        'enable'  , 'off',...
        'units'   , 'normalized',...
        'position', [0.8 0.02 0.15, 0.05],...
        'callback', @OK);
    
    % Adjust the settings of all UICONTROL's if the settings have been
    % changed before.
    set(best_of, 'string', BATCH.best_of);
    set(min_GAMs_popup, 'value', BATCH.min_GAMs);
    set(max_GAMs_popup, 'value', BATCH.max_GAMs);
    set(selected_bodies_list, 'string', BATCH.bodies);        
    % adjust all neccesary settings
    modify_controls('init');

    % NOW show window
    set(batch_window, 'visible', 'on')    
    % pause all other output
    uiwait(batch_window);

    %% All callback functions
    
    % abbreviate body names
    function abbreviations = abbreviate_names(which_ones, letters)
        % defaults
        if (nargin == 0), which_ones = 1:numel(body_names); letters = 1; end
        % form abbreviations
        abbreviations = cellfun(@(x)x(1:letters), body_names(which_ones), 'uniformoutput', false);
        % check for doubles
        for ii = 1:numel(abbreviations)
            % compare current string to all strings
            equals = strcmpi(abbreviations{ii}, abbreviations);
            % call this function recursively if some are equal
            if (nnz(equals) > 1)
                abbreviations(equals) = abbreviate_names(find(equals), letters+1);
            end
        end
    end
    
    % add a GAM-body to the list
    function add_body(varargin)
        % get the currently selected entry
        current_entry = get(available_bodies_list, 'value');
        % add the string to the list of selected bodies
        BATCH.bodies = get(selected_bodies_list, 'string');
        BATCH.bodies{end+1} = body_names{current_entry};
        set(selected_bodies_list, ...
            'string', BATCH.bodies,...
            'value' , length(BATCH.bodies));
        % enable or disable "OK", "Clear" and "Remove" buttons
        modify_controls;
    end % add body
    
    % remove a GAM-body from the list
    function remove_body(varargin)
        % get the currently selected body
        current_selection = get(selected_bodies_list, 'value');
        % remove the entry
        BATCH.bodies(current_selection) = [];
        set(selected_bodies_list, ...
            'string', BATCH.bodies,...
            'value' , current_selection-1);        
        % enable or disable "OK", "Clear" and "Remove" buttons
        modify_controls;
    end % remove body
    
    % clear the entire list
    function clear_all_bodies(varargin)
        % clear the strings
         set(selected_bodies_list, 'string', {});
         BATCH.bodies = {};
         % and modify all relevant controls
         modify_controls;
    end % clear_all_bodies
    
    % enable or disable appropriate controls
    function modify_controls(varargin)
        % get current stringlist
        list_of_selected_bodies = get(selected_bodies_list, 'string');
        % number of bodies in the list
        num_bodies = length(list_of_selected_bodies);        
        % enable or disable appropriate controls
        if (~isempty(list_of_selected_bodies))
            % first include all combinations
            if (nargin == 0)
                BATCH.combinations = all_combinations;
                BATCH.exclusions = false;
            end
            % modify all controls
            set(remove_button    , 'enable', 'on');
            set(OK_button        , 'enable', 'on');
            set(clear_button     , 'enable', 'on');
            set(combinations_text, 'enable', 'on');
            set(exclusion_list_button, 'enable', 'on');
            set(max_GAMs_popup, ...
                'enable', 'on',...
                'string', 1:num_bodies);           
            set(min_GAMs_popup, ...
                'enable', 'on',...
                'string', 1:num_bodies);            
            % make sure the selected values stay in range
            set(max_GAMs_popup, 'value', min(get(max_GAMs_popup, 'value'), num_bodies));
            set(min_GAMs_popup, 'value', min(get(min_GAMs_popup, 'value'), num_bodies));            
            % also show number of combinations
            show_number_of_optimizations;
        else
            % first include all combinations
            BATCH.combinations = {};
            % modify controls
            set(max_GAMs_popup,...
                'enable', 'off',...
                'string', 0,...
                'value' , 1);
            set(min_GAMs_popup,...
                'enable', 'off',...
                'string', 0,...
                'value' , 1);
            set(remove_button    , 'enable', 'off');
            set(OK_button        , 'enable', 'off');
            set(clear_button     , 'enable', 'off');
            set(combinations_text, 'enable', 'off');
            set(exclusion_list_button, 'enable', 'off');
            % also show number of combinations
            set(combinations_text(2), 'string', 0);
        end
        % toggle or untoggle "Exclude Combinations" button
        if BATCH.exclusions
            set(exclusion_list_button, 'value', 1);
        else
            set(exclusion_list_button, 'value', 0);
        end
    end % modify controls
    
    % calculate number of optimizations to perform
    function num_combinations = show_number_of_optimizations(varargin)
        % get current stringlist
        BATCH.bodies = get(selected_bodies_list, 'string');
        % calculate the number of possible combinations            
        BATCH.max_GAMs = get(max_GAMs_popup, 'value');
        BATCH.min_GAMs = get(min_GAMs_popup, 'value');
        % they might be swapped
        if (BATCH.min_GAMs > BATCH.max_GAMs)
            tmp = BATCH.max_GAMs; 
            BATCH.max_GAMs = BATCH.min_GAMs; 
            BATCH.min_GAMs = tmp; 
        end
        % calculate number of combinations
        if isempty(BATCH.combinations) || (nargin > 0 && ...
                any(varargin{1} == [min_GAMs_popup, max_GAMs_popup]))
            BATCH.combinations = all_combinations;
        end
        num_combinations = numel(BATCH.combinations);        
        BATCH.num_combinations = num_combinations;
        % multiply by "best of"
        BATCH.best_of = str2double(get(best_of, 'string'));
        % set new string        
        set(combinations_text(2), 'string', BATCH.best_of*num_combinations);
    end
    
    % OK-button callback
    function OK(varargin)
        % save settings
        settings.BATCH = BATCH;
        % set output argument
        accept = true;
        % resume MainWin
        uiresume(batch_window);             
        % and kill window
        delete(batch_window);        
    end % OK 
    
    % cancel button
    function cancel_button(varargin)
        % set output argument
        varargout{1} = false;
        % resume MainWin
        uiresume(batch_window);        
        % and kill window
        delete(batch_window);
    end % cancel_button
    
    % make all possible combinations
    function combinations = all_combinations(varargin)
        % initialize        
        combinations = {};
        % find indices of the selected bodies
        indices = zeros(1,numel(BATCH.bodies));
        for ii = 1:numel(BATCH.bodies)
            indices(ii) = find(strcmpi(body_names, BATCH.bodies{ii}), 1);
        end
        % all permutations
        all_combos = perms(indices);        
        % loop through the min/max GAMs 
        index = 0;
        for ii = BATCH.min_GAMs:BATCH.max_GAMs
            % unique combinations with these limits
            range_combos = unique(all_combos(:,1:ii),'rows');
            % convert to cell-array
            for jj = 1:size(range_combos,1)
                index = index + 1;
                combinations{index,1} = range_combos(jj,:);%#ok
            end
        end
    end % all combinations 
    
    %% Exclusion list window
    
    % sub-GUI of BATCH-GUI
    function exclusion_list(varargin)
           
        %% initialize
        
        % set a maximum of 800 unique combinations
        if (BATCH.num_combinations > 800)
            % show warning
            warndlg({'There are too many possible combinations to display in the exclusion list.'
                'Please consider reducing the number of swingby bodies or the number of swingby''s.'},...
                'Too many combinations');
            % and return
            return;
        end
        
        % initialize all combinations
        if ~BATCH.exclusions % if we've never been here before       
            BATCH.combinations = all_combinations;
        end
        
        % copy original
        modified_combinations = BATCH.combinations;
        
        %% create sub-GUI
                
        % make new GUI
        excl_list = figure(...
            'position', [scx scy 600 500],...     % should be fine
            'renderer', 'zbuffer', ...            % safest & fastest choice
            'name'    , 'Exclusion list for BATCH-optimization',... % "Figure 2" is just corny
            'numbertitle', 'off',...
            'visible' , 'on',...                  % don't show until the window is fully drawn
            'units'   , 'normalized',...          % better for resizing etc.
            'resize'  , 'on', ...                 % but just keep it resizable
            'menubar' , 'none', ...               % menubar is not needed
            'toolbar' , 'none', ...               % no toolbar
            'windowstyle' , 'modal',...           % keep it on top
            'DockControls', 'off',...             % force it to be non-dockable
            'color'   , environment.colors.window_bgcolor,... % use system-default colorscheme
            'defaultuicontrolfontsize', 8);       % force font size
        
        % static text
        uicontrol(...
            'units'   , 'normalized',...
            'style'   , 'text',...
            'position',[0.01 0.87 0.98, 0.10],...
            'string'  , 'List of all combinations to try')
        
        % build table
        combinations_table = uitable(...
            'units'   , 'normalized',...
            'position', [0.01, 0.4, 0.98, 0.52],...
            'columnname', {},...
            'rowname' , {},...
            'fontsize', 10,...
            'fontweight', 'bold',...
            'columnformat', {'char','char','char','char'},...
            'columnwidth' , {145, 145, 145, 145},...
            'CellSelectionCallback', @CellSelectionCallback);
        
        % legend of abbreviations
        % NOTE - this is just a simple column-printer. A more 
        % sophisticated one is preferable; one where the spacing 
        % is constant etc. But...too complicated for now :)
        legend_string = {}; spaces = repmat(' ', 1,20);
        for i = 1:2:numel(abbrv_names)-1
            legend_string = [legend_string; 
                {[abbrv_names{i+0}, ' = ', body_names{i+0}, spaces, ...
                  abbrv_names{i+1}, ' = ', body_names{i+1} ]}];%#ok
        end
        legend_panel = uipanel(...
            'units'   , 'normalized',...
            'position', [0.20 0.10 0.57 0.3],...
            'title'   , 'Legend');
            uicontrol(...
                'style'   , 'text',...
                'units'   , 'normalized',...
                'horizontalalignment', 'left',...
                'parent'  , legend_panel,...
                'position', [0.05 0.05 .90 .90],...
                'string'  , legend_string);
        
        % "REMOVE" button
        uicontrol(...
            'style'   , 'pushbutton',...
            'string'  , 'REMOVE >>',...
            'units'   , 'normalized',...
            'position', [0.8 0.3 0.15, 0.05],...
            'callback', @remove);
        
        % "Reset Table" button
        uicontrol(...
            'style'   , 'pushbutton',...
            'string'  , 'Reset table',...
            'units'   , 'normalized',...
            'position', [0.02 0.3 0.15, 0.05],...
            'callback', @reset_table);
        
        % "Cancel" button
        uicontrol(...
            'style'   , 'pushbutton',...
            'string'  , 'Cancel',...
            'units'   , 'normalized',...
            'position', [0.35 0.02 0.15, 0.05],...
            'callback', @Cancel);
        
        % "OK" button
        uicontrol(...
            'style'   , 'pushbutton',...
            'string'  , 'OK',...
            'units'   , 'normalized',...
            'position', [0.55 0.02 0.15, 0.05],...
            'callback', @OK);
        
        % Set initial table data
        set_table_data(modified_combinations);   
        selected_cells = [];
        % show GUI
        set(excl_list, 'visible', 'on');
        
        %% callbacks
        
        % NOTE: Using an external [modified_combinations] array is 
        % necessary to allow "Cancel" to work as expected; if you 
        % directly insert data in the BATCH-structure, you loose 
        % the original information, so Cancel won't work correctly.
        
        % change selected cell(s)
        function CellSelectionCallback(varargin)
            selected_cells = varargin{2}.Indices;
        end
        
        % "exclude" button callback
        function remove(varargin)  
            % check if some were actually selected
            if isempty(selected_cells), return, end
            % convert the cell-indices to indices in
            % [modified_combinations]
            indices = 4*(selected_cells(:,1)-1) + selected_cells(:,2);
            % indices might include empty cells
            indices(indices > numel(modified_combinations)) = [];
            % remove selected combination from the list
            modified_combinations(indices) = [];  
            % update table
            set_table_data(modified_combinations);
            % now we have exlusions
            BATCH.exclusions = true;
        end
        
        % adjust table
        function set_table_data(data)
            % first reshape data into proper format (4 columns)
            dat = {};            
            for ii = 1:4:size(data,1)                
                col1 = strcat(abbrv_names{data{ii+0}});
                col2 = ' '; col3 = ' '; col4 = ' ';
                if (ii + 1) <= size(data,1), col2 = strcat(abbrv_names{data{ii+1}}); end
                if (ii + 2) <= size(data,1), col3 = strcat(abbrv_names{data{ii+2}}); end
                if (ii + 3) <= size(data,1), col4 = strcat(abbrv_names{data{ii+3}}); end
                dat = [dat; {col1, col2, col3, col4}];%#ok 
            end
            % THEN update table
            set(combinations_table, 'data', dat);
        end
        
        % "Cancel" = don't save, kill window
        function Cancel(varargin)
            delete(excl_list);
        end
        
        % "OK" = first save, kill window
        function OK(varargin)
            % copy data
            BATCH.combinations     = modified_combinations;
            BATCH.num_combinations = size(BATCH.combinations, 1);
            % set new amount of optimizations to perform
            show_number_of_optimizations;
            % toggle or untoggle "Exclude combinations" button
            if BATCH.exclusions
                set(exclusion_list_button, 'value', 1);
            else
                set(exclusion_list_button, 'value', 0);
            end
            % and kill window
            delete(excl_list);
        end
        
        % reset table
        function reset_table(varargin)
            % reset structure
            BATCH.exclusions   = false; 
            BATCH.combinations = all_combinations;
            % reset table and copy data
            set_table_data(BATCH.combinations);
            modified_combinations = BATCH.combinations;
        end % reset table
        
    end % exclusion list
    
end % BATCH optimiziation
