function dirPath = thisMfileDir()
%% thisMfileDir
% Author: Erik Roberts
%
% Purpose: Gets current m-file dir
%
% Usage: dirPath = thisMfileDir()

stack = dbstack;
% dirPath = evalin('caller', 'fileparts(which(mfilename(''fullpath'')))');
dirPath = fileparts(which(stack(end).name));

end