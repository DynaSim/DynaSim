function [variables,pop_names]=SelectVariables(labels,var_strings)
% variables=SelectVariables(labels,var_strings)
% purpose: determine what variables to plot
% inputs:
%   labels - cell array of variable names
%   var_strings - string or cell array of strings specifying variables to plot
% outputs:
%   all labels matching specifications in var_strings
% examples:
%   labels={'pop1_v','pop1_iNa_m','pop1_iNa_h','pop2_v','pop2_av','time'};
%   var_strings=[];
%   var_strings='v';
%   var_strings='pop1';
%   var_strings='*_v';
%   var_strings='pop1_v';
%   var_strings='pop1_*';
%   var_strings='pop2_*';
if nargin<2
  var_strings=[];
end
if isempty(var_strings)
  % set default: all pops with state variable of first element of labels
  var=regexp(labels{1},'_.*$','match');
  % add wildcard
  if isempty(var)
    var_strings={'*'};
  else
    var_strings={['*' var{1}]};
  end
elseif ~iscell(var_strings)
  var_strings={var_strings};
end
% loop over cell array of variable indicators
variables={};
for i=1:length(var_strings)
  varstr=var_strings{i};
  % convert state variable into reg string to get variable for all pops
  if ~any(varstr=='*') && any(~cellfun(@isempty,regexp(labels,['_' varstr '$'])))
    varstr=['*_' varstr];
  end
  % convert any population name into reg string to get all variables in pop
  if ~any(varstr=='*') && any(~cellfun(@isempty,regexp(labels,['^' varstr '_'])))
    varstr=[varstr '_*'];
  end
  % add period to get all matches
  varstr=strrep(varstr,'*','.*');
  % find all matches
  matches=regexp(labels,['^' varstr '$'],'match');
  variables=cat(2,variables,matches{:});
end
if nargout>1
  pop_names={};
  for i=1:length(variables)
    name=regexp(variables{i},'^([a-zA-Z0-9]+)_','tokens','once');
    pop_names{end+1}=name{1};
  end
end
