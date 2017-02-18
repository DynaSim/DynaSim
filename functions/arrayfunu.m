

function varargout = arrayfunu(varargin)

    varargout = cell(1,nargout);
    [varargout{:}] = arrayfun(varargin{:},'UniformOutput',0);

end