classdef test_checkModel < matlab.unittest.TestCase
  properties
    unitTestDataPath = ds.getConfig('ds_unitTestData_path');
  end
  
  methods (Test,  TestTags = {'core'})
    function testCellIn(testCase)
      argin = {{'s=10; r=27; b=2.666';'dx/dt=s*(y-x)';'dy/dt=r*x-y-x*z';'dz/dt=-b*z+x*y'}};
      expectedOut = load(fullfile(testCase.unitTestDataPath,'checkModel_testCellIn'));
      expectedOut = expectedOut.model;
      
      testCase.verifyEqual(expectedOut, ds.checkModel(argin{:}) );
    end
  end
  
end