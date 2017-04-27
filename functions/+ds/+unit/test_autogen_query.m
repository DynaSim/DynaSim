classdef test_autogen_query < matlab.unittest.TestCase
  properties
    testDataPath = ds.getConfig('ds_testData_path');
  end
  
  properties (TestParameter)
    dataFileName = ds.unit.getAutogenFiles(false, true);
  end
  
  methods (Test,  TestTags = {'autogen', 'query'})
    function testCellIn(testCase, dataFileName)
      args = load(fullfile(testCase.testDataPath, 'autogen', dataFileName));
      expectedOut = args.argout;
      
      [~,filename] = fileparts2(dataFileName);
      fnName = strsplit(filename,'_autogen_');
      
      % TODO: handle localfunctions from query
      
      fh = str2func(fnName{1});
      
      [testOut{1:length(expectedOut)}] = feval(fh, args.argin{:});
      
      for ind = 1:length(expectedOut)
        testCase.verifyEqual(testOut{ind}, expectedOut{ind});
      end
      
    end
  end
  
end