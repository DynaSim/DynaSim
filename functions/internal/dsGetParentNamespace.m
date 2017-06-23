function parent = dsGetParentNamespace(namespace, varargin)
%GETPARENTNAMESPACE - determine parent namespace from namespace specified in namespace
%
% Usage:
%   parent = dsGetParentNamespace(namespace)
%
% Input:
%   - namespace: current namespace of object
%
% Output:
%   - parent: parent namespace containing the current namespace
%
% Examples:
%   parent=dsGetParentNamespace('pop')
%   parent=dsGetParentNamespace('pop_mech')
%   parent=dsGetParentNamespace('pop_pop')
%   parent=dsGetParentNamespace('pop_pop_mech')
%   parent=dsGetParentNamespace('mech')
%   parent=dsGetParentNamespace('')
%
% See also: dsPropagateNamespaces

%% auto_gen_test_data_flag argin
options = dsCheckOptions(varargin,{'auto_gen_test_data_flag',0,{0,1}},false);
if options.auto_gen_test_data_flag
  varargs = varargin;
  varargs{find(strcmp(varargs, 'auto_gen_test_data_flag'))+1} = 0;
  varargs(end+1:end+2) = {'unit_test_flag',1};
  argin = [{namespace}, varargs]; % specific to this function
end


if isempty(namespace) && isnumeric(namespace)
  namespace='';
end
if ~isempty(namespace) && namespace(end)=='_'
  namespace=namespace(1:end-1);
end
if ~isempty(namespace)
  parts=regexp(namespace,'_','split');
else
  parts=[];
end

switch length(parts)
  case 0                          % ''
    parent='global';
  case 1                          % pop or mech
    parent='';
  case 2
    if isequal(parts{1},parts{2}) % pop_pop
      parent='global';
    else                          % pop_mech
      parent=[parts{1} '_'];
    end
  case 3                          % pop_pop_mech
    parent=[parts{1} '_' parts{2} '_'];
  otherwise                       % a_b_c_d_...
    parent='';
    for i=1:length(parts)-1
      parent=[parent parts{i} '_'];
    end
end

%% auto_gen_test_data_flag argout
if options.auto_gen_test_data_flag
  argout = {parent}; % specific to this function
  
  dsUnitSaveAutoGenTestData(argin, argout);
end
