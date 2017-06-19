classdef dsUnitTest_autogen_all < matlab.unittest.TestCase
  properties
    unitTestDataPath = dsGetConfig('ds_unitTestData_path');
  end
  
  properties (TestParameter)
    dataFileName = dsUnitGetAutogenFiles();
  end
  
  methods (Test,  TestTags = {'autogen'})
    function testCellIn(testCase, dataFileName)
      args = load(fullfile(testCase.unitTestDataPath, 'autogen', dataFileName));
      expectedOut = args.argout;
      
      [~,filename] = fileparts2(dataFileName);
      fnName = strsplit(filename,'_autogen_');
      fh = str2func(fnName{1});
      
      [testOut{1:length(expectedOut)}] = feval(fh, args.argin{:});
      
      for ind = 1:length(expectedOut)
        testCase.verifyEqual(testOut{ind}, expectedOut{ind});
      end
      
    end
  end
  
end
