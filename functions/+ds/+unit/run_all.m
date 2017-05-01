%% Run all tests
% excluding query tests

import matlab.unittest.TestSuite
import matlab.unittest.selectors.HasTag
import matlab.unittest.plugins.CodeCoveragePlugin

%% Make Test Suite
fullSuite = TestSuite.fromPackage('ds.unit');
fullSuite = fullSuite.selectIf(~HasTag('query'));

%% code coverage runner
runner = TestRunner.withTextOutput;
runner.addPlugin(CodeCoveragePlugin.forFolder(fullfile(ds.getConfig('ds_root_path'), 'functions')))
runner.addPlugin(CodeCoveragePlugin.forPackage('ds'))

%% Run Test Suite
% result = run(fullSuite); % test suite without code covereage or parallel
result = runner.run(fullSuite); % runner for code coverage

%% Run Test Suite in Parallel
% runner = matlab.unittest.TestRunner.withTextOutput;
% result = runInParallel(runner,fullSuite); % runner in parallel, no code covereage