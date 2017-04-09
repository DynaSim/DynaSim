classdef test_autogen_all < matlab.unittest.TestCase
  properties
    testDataPath = ds.getConfig('ds_testData_path');
  end
  
  properties (TestParameter)
    dataFileName = lscell(fullfile(ds.getConfig('ds_testData_path'), 'autogen', '*_autogen_*'));
  end
  
  methods (Test,  TestTags = {'autogen'})
    function testCellIn(testCase, dataFileName)
      args = load(fullfile(testCase.testDataPath, 'autogen', dataFileName));
      expectedOut = args.argout;
      
      [~,filename] = fileparts2(dataFileName);
      fnName = strsplit(filename,'_autogen_');
      fh = str2func(fnName{1});
      
%       if ~isfunction(fh)
%         fh = str2func(['ds.' fnName{1}]);
%       end
      
      [actualOut{1:length(expectedOut)}] = feval(fh, args.argin{:});
      
      for ind = 1:length(expectedOut)
        testCase.verifyEqual(expectedOut{ind}, actualOut{ind});
      end
      
    end
  end
  
end