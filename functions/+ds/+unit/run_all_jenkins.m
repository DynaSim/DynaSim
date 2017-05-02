%% Run all tests
% excluding query tests

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.selectors.HasTag
import matlab.unittest.plugins.CodeCoveragePlugin
import edu.stanford.covert.test.Coverage

%% workspace
fprintf('Running from dir:%s\n',pwd);
[~,ws] = system('echo $WORKSPACE');
ws = strtrim(ws);
fprintf('Workspace:%s\n',ws);

%% Make Default Config
ds.makeDefaultConfigJenkins;

%% Make Test Suite
fullSuite = TestSuite.fromPackage('ds.unit');
fullSuite = fullSuite.selectIf(~HasTag('query'));

%% code coverage runner
% runner = TestRunner.withTextOutput;
% runner.addPlugin(CodeCoveragePlugin.forFolder(fullfile(ds.getConfig('ds_root_path'), 'functions')))
% runner.addPlugin(CodeCoveragePlugin.forPackage('ds'))

%% Run Test Suite
% results = run(fullSuite); % test suite without code covereage or parallel
% results = runner.run(fullSuite); % runner for code coverage

%% Run Test Suite in Parallel
runner = matlab.unittest.TestRunner.withTextOutput;
try
  results = runInParallel(runner,fullSuite); % runner in parallel, no code coverage
  display(results);

  %% XML Coverage Output
  testCoverageDir = fullfile(ws, 'testCoverage');
  mkdirSilent(testCoverageDir)
  reportPath = fullfile(testCoverageDir, 'dsAllTestCoverageJenkins.xml');
  report = Coverage( fullfile(ws, 'functions'), ws );
  report.exportXML(reportPath);
catch e
  disp(getReport(e,'extended'))
  exit(1);
end

exit(any([results.Failed]));