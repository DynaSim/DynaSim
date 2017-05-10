function output = cellfunu(varargin)
%CELLFUNU - Wraps cellfun, with uniform output set to zero
%
% Usage:
%   out = cellfunu (varargin)

output = cellfun(varargin{:},'UniformOutput',0);

end