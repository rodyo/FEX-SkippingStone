function varargout = progressOrbit(varargin)
% Little wrapper function to call the MEX function when it exists, 
% otherwise, fall back to the slower M-version
    
    if exist('progress_orbit', 'file')~=3         
        [varargout{1:nargout}] = progress_orbitM(varargin{:});
    else
        [varargout{1:nargout}] = progress_orbit(varargin{:});
    end
    
end
