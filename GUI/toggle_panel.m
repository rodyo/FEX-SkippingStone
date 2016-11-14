% Toggle a panel's children's enable states.
function states = toggle_panel(panels_in, action)

    for ii = 1:numel(panels_in)

        panel = panels_in(ii);

        switch lower(action)

            % Disable: get all current enable states of all children, and
            % store this in the Panel's UserData field
            case 'disable'

                % If there already is data in the field, it's probably a sub-panel
                % already managed separately; leave it be
                if ~isempty(get(panel, 'UserData'))
                    states = {};
                    return;
                end

                % Get all children's enable states
                states = {[], {}};

                children = get(panel   , 'children');
                types    = get(children, 'type');
                panelind = strcmp(types, 'uipanel') | strcmp(types, 'uibuttongroup');
                panels   = children( panelind);
                other    = children(~panelind);

                % Recurse for nested panels when needed
                for jj = 1:numel(panels)
                    new_states = toggle_panel(panels(jj), action);
                    if ~isempty(new_states)
                        states{1} = [states{1}; new_states{1}];
                        states{2} = [states{2}; new_states{2}];
                    end
                end

                states{1} = [states{1}; other(:)];
                states{2} = [states{2}; get(other, 'enable')];

                % Save states in mother panel
                set(panel, 'UserData', states);

                % Now disable all controls
                set(other, 'enable', 'off');

            case 'reset'

                % Get the previous states
                states = get(panel, 'UserData');

                % Apply the previous states to the relevant controls
                for jj = 1:numel(states{1})
                    set(states{1}(jj),...
                        'enable'  , states{2}{jj});
                end

                % Clear the UserData field
                set(panel, 'UserData', []);
                states = {};
        end

    end

end
