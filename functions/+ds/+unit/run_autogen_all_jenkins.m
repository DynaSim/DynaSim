%% Run all autogen tests
% result_localfn = runtests('ds.unit.test_autogen_all_localfn');
% result = runtests('ds.unit.test_autogen_all');

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.selectors.HasTag
import matlab.unittest.plugins.CodeCoveragePlugin
import edu.stanford.covert.test.Coverage

fprintf('Running from dir:%s\n',pwd);

%% Make Default Config
ds.makeDefaultConfigJenkins;

%% Make Test Suite
fullSuite = TestSuite.fromPackage('ds.unit');
fullSuite = fullSuite.selectIf(~HasTag('query'));
fullSuite = fullSuite.selectIf(HasTag('autogen'));

%% code coverage runner
% runner = TestRunner.withTextOutput;
% runner.addPlugin(CodeCoveragePlugin.forFolder(fullfile(ds_root_path, 'functions')))
% runner.addPlugin(CodeCoveragePlugin.forPackage('ds'))

%% Run Test Suite
% result = run(fullSuite); % test suite without code covereage or parallel
% result = runner.run(fullSuite); % runner for code coverage

%% Run Test Suite in Parallel
runner = matlab.unittest.TestRunner.withTextOutput;
try
  results = runInParallel(runner,fullSuite); % runner in parallel, no code coverage
  display(results);

  %% XML Coverage Output
  [~,ws] = system('echo $WORKSPACE');
  testCoverageDir = fullfile(ws, 'testCoverage');
  mkdirSilent(testCoverageDir)
  reportPath = fullfile(testCoverageDir, 'dsAllAutogenTestCoverageJenkins.xml');
  report = Coverage( fullfile(ws, 'functions') );
  report.exportXML(reportPath);
catch e
  disp(getReport(e,'extended'))
  exit(1);
end

exit(any([results.Failed]));