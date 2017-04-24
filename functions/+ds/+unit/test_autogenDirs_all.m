classdef test_autogenDirs_all < matlab.unittest.TestCase
  properties
    testDataPath = ds.getConfig('ds_testData_path');
  end
  
  properties (TestParameter)
    dataDirName = ds.unit.getAutogenDirs();
  end
  
  methods (Test,  TestTags = {'autogen'})
    function testCellIn(testCase, dataDirName)
      dataDirPath = fullfile(testCase.testDataPath, 'autogenDirs', dataDirName);
      
      % Make Temp Folder
      import matlab.unittest.fixtures.TemporaryFolderFixture
      tempDirFixture = testCase.applyFixture(TemporaryFolderFixture('WithSuffix', ['_' dataDirName]));
      tempDirPath = tempDirFixture.Folder;
      
      args = load(fullfile(dataDirPath, 'args.mat'));
      expectedOut = args.argout;
      
      fnName = strsplit(dataDirName,'_autogen_');
      fh = str2func(fnName{1});
      
      % TODO: use tempDirPath for actual outputs
      
      [actualOut{1:length(expectedOut)}] = feval(fh, args.argin{:});
      
      % actualOut fig handles
      % TODO: generalize fig handle checking
      if all(isValidFigHandle(actualOut{1}))
        handles = actualOut{1};
        
        ds.unit.save_figHandles( handles, tempDirPath )
        
        close all
      end
      
      % test output args
      for ind = 1:length(expectedOut)
        testCase.verifyEqual(expectedOut{ind}, actualOut{ind});
      end
      
      % test output files
      outputDir = fullfile(dataDirPath, 'output');
      if exist(outputDir, 'dir')
        expectedOutputFiles = lscell(outputDir);
        actualOutputFiles = lscell(tempDirPath);
        
        % test same num files
        testCase.assertLength(length(actualOutputFiles), length(expectedOutputFiles));
        
        % test each file
        for iFile = 1:length(expectedOutputFiles)
          thisExpectedOutputFilePath = fullfile(outputDir, expectedOutputFiles{iFile});
          
          thisActualOutputFilePath = fullfile(tempDirPath, actualOutputFiles{iFile});
          [~,~,ext] = fileparts2(thisActualOutputFilePath);
          switch ext
            case '.fig'
              % compareFigFiles 
%               testCase.verifyTrue( ~any(cellfun(@isempty, compareFigFiles(actualOutputFiles{iFile}, expectedOutputFiles{iFile}))) );
              testCase.verifyEmpty( compareFigFiles(thisActualOutputFilePath, thisExpectedOutputFilePath, true) );
            case '.mat'
              % compare loaded structs
              thisActual = load(thisActualOutputFilePath);
              thisExpected = load(thisExpectedOutputFilePath);
              testCase.verifyEqual(thisActual, thisExpected);
            otherwise
              % do diff
              testCase.verifyFalse( logical(system(sprintf('diff %s %s', thisActualOutputFilePath, thisExpectedOutputFilePath))) );
          end
        end
      end
      
    end
  end
  
end