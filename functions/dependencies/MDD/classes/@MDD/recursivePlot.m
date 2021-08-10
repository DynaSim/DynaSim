function varargout = recursivePlot(varargin)
    % Alias for recursiveFunc
    
    [varargout{1:nargout}] = recursiveFunc(varargin{:});
    
end