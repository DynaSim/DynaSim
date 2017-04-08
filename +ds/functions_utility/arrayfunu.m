function varargout = arrayfunu(varargin)
%ARRAYFUNU - Wraps arrayfun, with uniform output set to zero
%
% Usage:
%   out = arrayfunu (varargin)

varargout = cell(1,nargout);
[varargout{:}] = arrayfun(varargin{:},'UniformOutput',0);

end
