function progress_bar(progress, string)
% Last changed 21/Jun/2009 (Rody)
   
% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com
% Affiliation: Delft University of Technology
%              Partly performed at ISAS/JAXA

    % initialize persistents & globals
    global MainWin % handle to main window
    persistent BGCOLOR HANDLE
    
    % define progress polygon and text
    if (nargin == 2)
        
        % define the polygon and line
        xpatch = [0, progress, progress, 0];
        ypatch = [0, 0, 1, 1]; 
        
        % add percentage to the string
        progress = max(0, min(1, progress));
        string = [string, ' (', num2str(round(progress*100)), '%)'];
        
        % define the color for the progress polygon        
        bgcolor = BGCOLOR + 0.1;
        bgcolor(bgcolor > 1) = 1;
    
    % (reset progress bar to defaults)
    elseif (nargin == 1)
        
        % first call is single cell array
        if ~isa(progress, 'char') 
            % set persistent variables
            BGCOLOR = progress{1};
            HANDLE  = progress{2};        
        end
        
        % subsequent calls are empty strings
        
        % define the polygon and line
        xpatch = [0, 1, 1, 0];
        ypatch = [0, 0, 1, 1];          
        % reset the string to "Idle"
        string = '(Idle)';
        
        % set backgournd color
        bgcolor = BGCOLOR;
        
    end
    
    % when this function is called when the Main Window is not running, do
    % nothing and simply return 
    if isempty (HANDLE) || isempty(BGCOLOR) || ~ishandle(HANDLE)
        return
    end
    
    % set the progress bar axes as the current axes and clear them
    set(0, 'currentfigure', MainWin);
    set(MainWin, 'currentAxes', HANDLE); cla;
    
    % lines to redraw the box 
    % (gets overwritten by the progress polygon)
    xline  = [1 0 0 1 1];
    yline  = [0 0 1 1 0];
    
    % draw the polygon, and re-draw the box
    patch(xpatch, ypatch, bgcolor, 'EdgeColor', bgcolor);
    l = line(xline,yline);   
    set(l, 'color', 'k')
    
    % set the text    
    text(0.5, 0.5, string, ...        
        'FontSize' , 9, ...
        'EraseMode', 'xor',...
        'color'    , 'k', ...
        'HorizontalAlignment', 'center');   
    
    % force the axes to be drawn (even in background)
    drawnow
    
end

