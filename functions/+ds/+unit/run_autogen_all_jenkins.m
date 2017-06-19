%% Run all autogen tests
% result_localfn = runtests('dsUnitTest_autogen_all_localfn');
% result = runtests('dsUnitTest_autogen_all');

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

%% fix paths
fprintf('Fixing paths.\n');
rmPathVar('ds');
addpath(genpath(ws));

%% Make Default Config
dsMakeDefaultConfigJenkins;

%% Rename autogen_newSave to autogen
finalDir = fullfile(dsGetConfig('ds_unitTestData_path'), 'autogen');
if ~exist(finalDir, 'dir')
  newDir = fullfile(dsGetConfig('ds_unitTestData_path'), 'autogen_newSave');
  movefile(newDir, finalDir);
end

%% Rename autogenDirs_newSave to autogenDirs
finalDir = fullfile(dsGetConfig('ds_unitTestData_path'), 'autogenDirs');
newDir = fullfile(dsGetConfig('ds_unitTestData_path'), 'autogenDirs_newSave');
if ~exist(finalDir, 'dir') && exist(newDir, 'dir')
  movefile(newDir, finalDir);
end

%% Make Test Suite
fullSuite = TestSuite.fromPackage('dsUnit');
fullSuite = fullSuite.selectIf(~HasTag('query'));
fullSuite = fullSuite.selectIf(HasTag('autogen'));

%% code coverage runner
% runner = TestRunner.withTextOutput;
% runner.addPlugin(CodeCoveragePlugin.forFolder(fullfile(ds_root_path, 'functions')))
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
  mkdirSilent(testCoverageDir);
  reportPath = fullfile(testCoverageDir, 'dsAllAutogenTestCoverageJenkins.xml');
  report = Coverage( fullfile(ws, 'functions'), ws );
  report.exportXML(reportPath);

  %% get coverage percentage
  coverPercent = report.stats.lineRate;
  if isnan(coverPercent)
    coverPercent = 0;
  end
  coverPercentStr = sprintf('%.f', coverPercent);
  fprintf('Test Coverage (%% of code lines): %s\n', coverPercentStr);

  %% make coverage svg
  filePath = '/home/erik/Dropbox/research/dsJenkinsBadge/coverage.svg';
  system(['python /home/erik/Dropbox/Programming/Python/universal_coverage_badge/__main__.py -percent ' coverPercentStr ' -o ' filePath ' -f True']);
catch e
  disp(getReport(e,'extended'))
  exit(1);
end

exit(any([results.Failed]));
