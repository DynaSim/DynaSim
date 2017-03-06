

% This function is now part of the main classdef file, nDDict.m


% This code is for enabling the use of string  data types as well.
% Disabling for now, since string only works in 2016b and backwards
% compatibility is complicated.
% function varargout = string(varargin)
%     % Overload string function for backwards compatibility
%     if exist('string')
%         [varargout{1:nargout}] = builtin('string',varargin{:});
%     else
%         c = varargin{1};
%         if iscell(c)
%             is_cellarray_of_chars(c);
%             out = c;
%         elseif ischar(c)
%             out = {c};
%         else
%             error('Unknown input type');
%         end
%         
%         varargout{1} = out;
%     end
% end
% 
% function varargout = isstring(varargin)
%     if exist('isstring','builtin')
%         [varargout{1:nargout}] = builtin('isstring',varargin{:});
%     else
%         c = varargin{1};
%         if iscell(c)
%             out = is_cellarray_of_chars(c);
%         elseif ischar(c)
%             out = true;
%         else
%             out=false;
%         end
%         varargout{1} = out;
%     end
% end
