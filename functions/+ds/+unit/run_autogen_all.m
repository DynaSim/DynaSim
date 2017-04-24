%% Run all autogen tests
% result_localfn = runtests('ds.unit.test_autogen_all_localfn');
% result = runtests('ds.unit.test_autogen_all');

import matlab.unittest.TestSuite
import matlab.unittest.selectors.HasTag

%% Make Test Suite
fullSuite = TestSuite.fromPackage('ds.unit');
fullSuite = fullSuite.selectIf(~HasTag('query'));
fullSuite = fullSuite.selectIf(HasTag('autogen'));

%% Run Test Suite
result = run(fullSuite);