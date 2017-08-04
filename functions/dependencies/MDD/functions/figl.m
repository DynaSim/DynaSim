function varargout = figl(varargin)
    % Create a large figure

    [varargout{1:nargout}] = figure('Units','normalized','Position',[0,0,1,1],varargin{:});

end