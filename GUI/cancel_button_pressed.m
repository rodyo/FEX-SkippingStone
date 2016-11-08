function yes_or_no = cancel_button_pressed(varargin)
%CANCEL_BUTTON_PRESSED              Returns the state of the cancel-button
%
% The function returns [true] if the Cancel-button has been pressed, and
% [false] if not.
%
% Note that you do not need to provide any arguments; it has been properly
% initialized by MAIN(). If the Main Window is NOT running, this function
% simply returns [false] (handy for debugging).

    % initialize persistent variable
    persistent handle

    % Its value is the handle to the cancel-button.
    % The handle is given upon initialization of the Main Window
    if (nargin > 0)
            handle = varargin{1};

    % This is convenient, because subsequent calls do not need to
    % know the handle. If called from some other function, its state
    % can directly be returned:
    else
        % if the Main Window is NOT running, return false
        if isempty(handle) || ~ishandle(handle)
            yes_or_no = false;
        % if it IS running, simply get the Cancel-button's state
        else
            yes_or_no = get(handle, 'UserData');
        end
    end
end
