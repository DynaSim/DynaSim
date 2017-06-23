%% Run non-dir non-local autogen tests

%% Rename autogen_newSave to autogen
finalDir = fullfile(ds.getConfig('ds_unitTestData_path'), 'autogen');
if ~exist(finalDir, 'dir')
  newDir = fullfile(ds.getConfig('ds_unitTestData_path'), 'autogen_newSave');
  movefile(newDir, finalDir);
end

%% Run tests
result = runtests('dsUnit.test_autogen_all', 'UseParallel', true);