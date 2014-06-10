function varargout = example_postprocessor(varargin)
% EXAMPLE
    
    %% return info & function handles    
        
    % build the infostructure
    pp_info.name = 'Example plugin';
    pp_info.GUI_function_handle = @(value)GUI_interaction(value);
    pp_info.plot_function_handle = @(varargin) plot_fcn(varargin{:});
    pp_info.function_handle = @(obj, x0, fval, LB, UB)...
                                pp_fcn(obj, x0, fval, LB, UB);    
    % return the info & return
    if (nargin == 0), varargout{1} = pp_info; return; end
        
    %% What should we do when it's selected in the GUI? 
    
    function GUI_interaction(value)%#ok
        % initalize        
        global MainWin arrival_tab
        handles = getappdata(MainWin, 'handles');
        % explanatory text
        handles.tab(arrival_tab).post_processing.additional_controls.panel(1) = uipanel(...
            'parent'    , handles.tab(arrival_tab).post_processing_panel,...
            'position'  , [0.01 0.2 0.98 0.6],...
            'bordertype', 'none');
        uicontrol(...
            'parent'  , handles.tab(arrival_tab).post_processing.additional_controls.panel(1),...
            'units'   , 'normalized',...
            'style'   , 'text',...
            'horizontalalignment', 'left',...
            'position', [0.01 0.77 0.98 0.23],...
            'string'  , {'This is a simple example of how to write a post-processor '
            'plugin; it doesn''t do any calculations.'});
        % finalize
        setappdata(MainWin, 'handles', handles);
    end % GUI function
    
    %% What should we do when any results are to be plotted? 
    
    function plot_fcn(varargin)
        % it doesn't do anything
    end % plot function
    
    %% What is the actual post-processing routine?
    
    %if (settings.postprocessing.post_processor == 2)
    function [postprocessor_data, solution, fval] = ...
              pp_fcn(obj, x0, fval, LB, UB)%#ok
        % it doesn't do anything; return "nothing"
        postprocessor_data = [];
        solution           = x0;
    end % pp_fcn
    
end % example post-processor plugin
