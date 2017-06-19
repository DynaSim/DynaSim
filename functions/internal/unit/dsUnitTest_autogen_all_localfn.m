classdef dsUnitTest_autogen_all_localfn < matlab.unittest.TestCase
  properties
    unitTestDataPath = dsGetConfig('ds_unitTestData_path');
  end
  
  properties (TestParameter)
    dataFileName = dsUnitGetAutogenFiles(true, false);
  end
  
  methods (Test,  TestTags = {'autogen'})
    function testCellIn(testCase, dataFileName)
      args = load(fullfile(testCase.unitTestDataPath, 'autogen', dataFileName));
      expectedOut = args.argout;
      
      % get fn and local fn names
      [~,filename] = fileparts2(dataFileName);
      dataNameSplit = strsplit(filename,'_autogen_');
      fnNameCat = dataNameSplit{1};
      splitName = strsplit(fnNameCat,'__');
      fnName = splitName{1};
      localfFnName = splitName{2};
      
      % Get local function handles from fn
      fnHandle = str2func(fnName);
      localFnHandles = feval(fnHandle);
      localFnStrs = cellfun(@func2str,localFnHandles, 'uni',false);
      
      % Get this local fn handle
      localFnHandle = localFnHandles{strcmp(localFnStrs, localfFnName)};
      
      [testOut{1:length(expectedOut)}] = feval(localFnHandle, args.argin{:});
      
      for ind = 1:length(expectedOut)
        testCase.verifyEqual(testOut{ind}, expectedOut{ind});
      end
      
    end
  end
  
end
