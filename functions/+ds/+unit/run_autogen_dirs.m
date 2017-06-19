%% Run autogen dirs tests

%% Rename autogenDirs_newSave to autogenDirs
finalDir = fullfile(dsGetConfig('ds_unitTestData_path'), 'autogenDirs');
if ~exist(finalDir, 'dir')
  newDir = fullfile(dsGetConfig('ds_unitTestData_path'), 'autogenDirs_newSave');
  movefile(newDir, finalDir);
end

%% Run tests
result = runtests('dsUnitTest_autogenDirs_all', 'UseParallel', true);