function dirList = lscell(arg, removePathBool, recursiveBool)
%% lscell
% Author: Erik Roberts
%
% Purpose: returns cell matrix of strings with results from call to ls/dir
%
% Usage: dirList = lscell()
%        dirList = lscell(arg)
%        dirList = lscell(arg, removePathBool)
%
% Inputs (optional):
%   arg: argument to ls
%   removePathBool: logical of whether to remove the path before the files/dirs

% args
if nargin < 2 || isempty(removePathBool)
  removePathBool = true; %defaults to true
end
if nargin < 3 || isempty(recursiveBool)
  recursiveBool = false; %defaults to false
end
if ~nargin || isempty(arg)
  if ~recursiveBool
    arg = {};
  else
    arg = {'.'};
  end
else
  arg = {arg};
end
  
if ~recursiveBool
  try
    dirList = ls(arg{:});
  catch % no arg files found
    dirList = {};
    return
  end
else % recursiveBool
  if (ismac || isunix)
    [~, dirList] = system(['find ' arg{:}]);
  elseif ispc
  end
end

if (ismac || isunix)
  dirList = strsplit(dirList);
  if isempty(dirList{end})
    dirList(end) = [];
  end
elseif ispc
  dirList = cellstr(dirList);
  if strcmp(dirList{1},'.')
    dirList(1:2) = []; %remove leading . and ..
  end
end

if isempty(dirList) || isempty(dirList{1})
  dirList = {};
end

if exist('removePathBool','var') && removePathBool
  dirList = removePath(dirList);
end
  
  % sub functions
  function dirList = removePath(dirList)
    for k = 1:length(dirList)
      thisDir = dirList{k};
      [~,name,ext] = fileparts(thisDir);
      dirList{k} = [name,ext];
    end
  end
end
