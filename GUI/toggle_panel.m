function states = toggle_panel(panel, action)
    
    switch lower(action)
        case 'disable'
            
            % Save the uicontrol states; easier to reenable that way
            states = {[], {}};
            
            children = get(panel   , 'children');
            types    = get(children, 'type');
            panelind = strcmp(types, 'uipanel');
            panels   = children( panelind);
            other    = children(~panelind);
            
            % Recurse for nested panels when needed
            for ii = 1:numel(panels)
                new_states = toggle_panel(panels(ii), action);
                states{1}  = [states{1}; new_states{1}];
                states{2}  = [states{2}; new_states{2}];
                % (no need for this on all levels)
                set(panels(ii), 'UserData', []);
            end
            
            states{1} = [states{1}; other(:)];
            states{2} = [states{2}; get(other, 'enable')];
            
            % Save states in mother panel
            set(panel, 'UserData', states);
            
            % Now disable all controls
            set(other, 'enable', 'off');
            
        case 'reset'
            
            states = get(panel, 'UserData'); 
            
            for ii = 1:numel(states{1})
                set(states{1}(ii),...
                    'enable'  , states{2}{ii});
            end
            
            states = {};
    end
    
end
