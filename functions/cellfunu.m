

function varargout = cellfunu (varargin)
    %out = cellfunu (varargin)
    % Cellfun with uniform output set to zero
    
    varargout = cell(1,nargout);
    [varargout{1:nargout}] = cellfun(varargin{:},'UniformOutput',0);

end