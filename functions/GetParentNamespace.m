function parent = GetParentNamespace(namespace)
%GETPARENTNAMESPACE - determine parent namespace from namespace specified in namespace
%
% Usage:
%   parent = GetParentNamespace(namespace)
%
% Input:
%   - namespace: current namespace of object
%
% Output:
%   - parent: parent namespace containing the current namespace
%
% Examples:
%   parent=GetParentNamespace('pop')
%   parent=GetParentNamespace('pop_mech')
%   parent=GetParentNamespace('pop_pop')
%   parent=GetParentNamespace('pop_pop_mech')
%   parent=GetParentNamespace('mech')
%   parent=GetParentNamespace('')
%
% See also: PropagateNamespaces

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
