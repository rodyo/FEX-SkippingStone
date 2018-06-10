function varargout = progressOrbit(varargin)
% Little wrapper function to call the MEX function when it exists, 
% otherwise, fall back to the slower M-version
    
    if exist('propagate_orbit', 'file')~=3         
        [varargout{1:nargout}] = propagate_orbitM(varargin{:});
    else
        [varargout{1:nargout}] = propagate_orbit(varargin{:});
    end
    
end
