%% Run non-dir non-local autogen tests

%% Rename autogen_newSave to autogen
finalDir = fullfile(dsGetConfig('ds_unitTestData_path'), 'autogen');
if ~exist(finalDir, 'dir')
  newDir = fullfile(dsGetConfig('ds_unitTestData_path'), 'autogen_newSave');
  movefile(newDir, finalDir);
end

%% Run tests
result = runtests('dsUnitTest_autogen_all', 'UseParallel', true);