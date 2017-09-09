classdef SkippingStoneGui < handle
    
    %% Properties
    % =====================================================================
    
    properties        
        handles        
    end
    
    properties (Constant, Access=private)
        
        mainwin_tag = 'SS_main_window';
    end
    
    
    %% Events
    % =====================================================================
    
    events
    end
    
    
    %% Methods
    % =====================================================================
    
    methods
    end
    
    methods (Access = private)
        
        % GUI is singleton 
        function obj = SkippingStoneGui(varargin)
            
        end
        
    end
    
    methods (Static)
        
        % Singleton constructor
        function obj = instance(varargin)
            
            persistent local_obj
            
            if isempty(local_obj)
                
                % The persistent may get cleared; isempty is not a
                % sufficient criterion for new-instance creation
                mw = findall(0, 'Tag', SkippingStoneGui.mainwin_tag);
                if isempty(mw)                
                    local_obj = SkippingStoneGui(varargin{:}); 
                else
                    local_obj = mw;
                end
            end            
            
            obj = local_obj;          
        end
        
    end
    
end 
