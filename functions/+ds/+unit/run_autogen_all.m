%% Run all autogen tests
% result_localfn = runtests('ds.unit.test_autogen_all_localfn');
% result = runtests('ds.unit.test_autogen_all');

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.selectors.HasTag
import matlab.unittest.plugins.CodeCoveragePlugin

%% Rename autogen_newSave to autogen
finalDir = fullfile(ds.getConfig('ds_unitTestData_path'), 'autogen');
if ~exist(finalDir, 'dir')
  newDir = fullfile(ds.getConfig('ds_unitTestData_path'), 'autogen_newSave');
  movefile(newDir, finalDir);
end

%% Rename autogenDirs_newSave to autogenDirs
finalDir = fullfile(ds.getConfig('ds_unitTestData_path'), 'autogenDirs');
newDir = fullfile(ds.getConfig('ds_unitTestData_path'), 'autogenDirs_newSave');
if ~exist(finalDir, 'dir') && exist(newDir, 'dir')
  movefile(newDir, finalDir);
end

%% Make Test Suite
fullSuite = TestSuite.fromPackage('ds.unit');
fullSuite = fullSuite.selectIf(~HasTag('query'));
fullSuite = fullSuite.selectIf(HasTag('autogen'));

%% code coverage runner
runner = TestRunner.withTextOutput;
runner.addPlugin(CodeCoveragePlugin.forPackage('ds'))

%% Run Test Suite
% result = run(fullSuite); % test suite without code covereage or parallel
result = runner.run(fullSuite); % runner for code coverage

%% Run Test Suite in Parallel
% runner = matlab.unittest.TestRunner.withTextOutput;
% result = runInParallel(runner,fullSuite); % runner in parallel, no code covereage