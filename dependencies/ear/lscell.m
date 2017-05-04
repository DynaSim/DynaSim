function dirList = lscell(arg, removePathBool)
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

if ~exist('removePathBool', 'var')
  removePathBool = true; %defaults to true
end

if nargin == 0
  dirList = ls;
else
  dirList = ls(arg);
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

if isempty(dirList{1})
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
