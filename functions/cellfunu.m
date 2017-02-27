function varargout = cellfunu (varargin)
%CELLFUNU - Wraps cellfun, with uniform output set to zero
%
% Usage:
%   out = cellfunu (varargin)

varargout = cell(1,nargout);
[varargout{1:nargout}] = cellfun(varargin{:},'UniformOutput',0);

end
