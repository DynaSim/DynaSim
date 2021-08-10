function dsRemoveLockFiles(dirPath)
% Purpose: Remove hidden lock files that can be left during errors and prevent
% certain functions from proceeding
%
% Author: Erik Roberts

if ~nargin
  dirPath = '';
end

removePathBool = false;
lockFiles = lscell( fullfile(dirPath, '.lock*'), removePathBool);

if ~isempty(lockFiles)
  cellfun(@delete, lockFiles);
end

end