function exist_mkdir(directoryPath)
%% exist_mkdir
% Author: Erik Roberts
%
% Purpose: Makes dir if doesn't exist
%
% Usage: exist_mkdir(directoryPath)

if ~exist(directoryPath, 'dir')
  mkdir(directoryPath);
end

end