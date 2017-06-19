%% Run all custom tests
% i.e. all tests that aren't autogen
import matlab.unittest.TestSuite
import matlab.unittest.selectors.HasTag

%% Make Test Suite
fullSuite = TestSuite.fromPackage('dsUnit');
fullSuite = fullSuite.selectIf(~HasTag('query'));
fullSuite = fullSuite.selectIf(~HasTag('autogen'));

%% Run Test Suite
result = run(fullSuite);