%% Test Class Definition
% References:
%  https://www.mathworks.com/help/matlab/matlab_prog/author-class-based-unit-tests-in-matlab.html
%  https://www.mathworks.com/help/matlab/matlab_prog/write-setup-and-teardown-code-using-classes.html
%  https://www.mathworks.com/help/matlab/matlab_prog/write-test-using-shared-fixtures.html
%  https://www.mathworks.com/help/matlab/matlab_prog/create-basic-parameterized-test.html

classdef template < matlab.unittest.TestCase
  
  %   properties % shared properties
  %     variable
  %   end
  %
  %   properties (TestParameter) % parametric cases
  %     type = {'single','double','uint16'};
  %     level = struct('small', 2,'medium', 4, 'large', 6);
  %     side = struct('small', 9, 'medium', 81,'large', 729);
  %   end
  %
  %   methods (TestMethodSetup)
  %     function createFixture(testCase)
  %       testCase.variable = figure;
  %       testCase.addTeardown(@close, testCase.variable)
  %     end
  %   end
  
  %% Test Method Block
%   methods (Test,  TestTags = {'Tag'})
%     function test1(testCase)
%       testCase.verifyEqual(testOut, expectedOut);
%     end
    
%     function parametricTest1(testCase, type, level)
%     end
%   end
end