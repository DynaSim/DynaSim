function keyval = Options2Keyval(options)
%OPTIONS2KEYVAL - Convert from options structure to a list of key/value pairs.
%
% Usage:
%   keyvals = Options2Keyval(options)
%
% Inputs:
%   - options: options structure to convert
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
