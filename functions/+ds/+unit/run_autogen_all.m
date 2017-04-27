%% Run all autogen tests
% result_localfn = runtests('ds.unit.test_autogen_all_localfn');
% result = runtests('ds.unit.test_autogen_all');

import matlab.unittest.TestSuite
import matlab.unittest.selectors.HasTag

%% Rename autogen_newSave to autogen
finalDir = fullfile(ds.getConfig('ds_testData_path'), 'autogen');
if ~exist(finalDir, 'dir')
  newDir = fullfile(ds.getConfig('ds_testData_path'), 'autogen_newSave');
  movefile(newDir, finalDir);
end

%% Rename autogenDirs_newSave to autogenDirs
finalDir = fullfile(ds.getConfig('ds_testData_path'), 'autogenDirs');
newDir = fullfile(ds.getConfig('ds_testData_path'), 'autogenDirs_newSave');
if ~exist(finalDir, 'dir') && exist(newDir, 'dir')
  movefile(newDir, finalDir);
end

%% Make Test Suite
fullSuite = TestSuite.fromPackage('ds.unit');
fullSuite = fullSuite.selectIf(~HasTag('query'));
fullSuite = fullSuite.selectIf(HasTag('autogen'));

%% Run Test Suite
result = run(fullSuite);