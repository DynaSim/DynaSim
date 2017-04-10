%% Run all tests
% excluding query tests

import matlab.unittest.TestSuite
import matlab.unittest.selectors.HasTag

%% Make Test Suite
fullSuite = TestSuite.fromPackage('ds.unit');
fullSuite = fullSuite.selectIf(~HasTag('query'));

%% Run Test Suite
result = run(fullSuite);