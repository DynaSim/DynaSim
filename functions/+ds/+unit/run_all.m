%% Run all tests
% excluding query tests

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.selectors.HasTag
import matlab.unittest.plugins.CodeCoveragePlugin
import edu.stanford.covert.test.Coverage

%% Make Test Suite
fullSuite = TestSuite.fromPackage('dsUnit');
fullSuite = fullSuite.selectIf(~HasTag('query'));

%% code coverage runner
runner = TestRunner.withTextOutput;
runner.addPlugin(CodeCoveragePlugin.forFolder(fullfile(dsGetConfig('ds_root_path'), 'functions')))
runner.addPlugin(CodeCoveragePlugin.forPackage('ds'))

%% Run Test Suite
% result = run(fullSuite); % test suite without code covereage or parallel
result = runner.run(fullSuite); % runner for code coverage

%% Run Test Suite in Parallel
% runner = matlab.unittest.TestRunner.withTextOutput;
% result = runInParallel(runner,fullSuite); % runner in parallel, no code coverage

%% XML Coverage Output
testCoverageDir = fullfile(dsGetConfig('ds_root_path'), 'testCoverage');
mkdirSilent(testCoverageDir)
reportPath = fullfile(testCoverageDir, 'dsAllTestCoverage.xml');
report = Coverage( fullfile(dsGetConfig('ds_root_path'), 'functions') );
report.exportXML(reportPath);