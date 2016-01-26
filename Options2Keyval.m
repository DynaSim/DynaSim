function keyval = Options2Keyval(options)
% keyvals = Options2Keyval(options)
% Purpose:
%   Convert from options structure to a list of key/value pairs.
%
% Parameters:
%   parms - options structure to convert
%
% See also: CheckOptions

% Grab the field names
fields = fieldnames(options);
keyval   = {};

% Loop over the field names, grab the value, and append both to the output
% list of key/value pairs
for i=1:length(fields)
  keyval{end+1} = fields{i};
  keyval{end+1} = options.(fields{i});
end;
