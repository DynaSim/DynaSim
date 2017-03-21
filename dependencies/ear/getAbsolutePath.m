function pathOut = getAbsolutePath(pathIn)
%% getAbsolutePath
%
% Purpose: Get absolute path from given absolute or relative path
%
% Usage: pathOut = getAbsolutePath(pathIn)
%
% Author: Erik Roberts

if ((isunix || ismac) && pathIn(1) == filesep) || (ispc && pathIn(2) == ':')
  % absolute
  pathOut = pathIn;
elseif pathIn(1) == '~'
  % relative to home
  [~,homePath] = system('echo $HOME');
  pathOut = fullfile(homePath(1:end-1), pathIn(2:end)); % need to remove return char
else
  % relative
  pathOut = fullfile(pwd, pathIn);
end

end