%% Run template
% References:
%   https://www.mathworks.com/help/matlab/matlab_prog/run-tests-for-various-workflows.html
%   https://www.mathworks.com/help/matlab/matlab_prog/create-simple-test-suites.html
%   https://www.mathworks.com/help/matlab/matlab_prog/tag-unit-tests.html

%% Imports
import matlab.unittest.TestSuite
import matlab.unittest.selectors.HasTag

%% Test Suites
% Examples:
%   fullSuite(1) = TestSuite.fromClass(?test_class);
%   fullSuite = TestSuite.fromFolder(pwd);
%   fullSuite = TestSuite.fromPackage('package');

%% Run Test Suite
% result = run(fullSuite);

%% Run Just Tagged Tests
% testNamesCell{fullSuite.Name}';
% results = runtests(testNamesCell,'Tag','Tag1');
%or
% taggedSuite = fullSuite.selectIf(HasTag('Tag1'));
% result = run(taggedSuite);