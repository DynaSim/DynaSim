classdef dsUnitTest_autogenDirs_all < matlab.unittest.TestCase
  properties
    unitTestDataPath = dsGetConfig('ds_unitTestData_path');
  end
  
  properties (TestParameter)
    dataDirName = dsUnitGetAutogenDirs();
  end
  
  methods (Test,  TestTags = {'autogen'})
    function testCellIn(testCase, dataDirName)
      dataDirPath = fullfile(testCase.unitTestDataPath, 'autogenDirs', dataDirName);
      
      % Make Temp Folder
      import matlab.unittest.fixtures.WorkingFolderFixture
      tempDirFixture = testCase.applyFixture(WorkingFolderFixture('WithSuffix', ['_' dataDirName]));
      tempDirPath = tempDirFixture.Folder;
%       tempOutputDir = fullfile(tempDirPath, 'output');
%       mkdirSilent(tempOutputDir);
      
      args = load(fullfile(dataDirPath, 'args.mat'));
      expectedOut = args.argout;
      
      fnName = strsplit(dataDirName,'_autogen_');
      fh = str2func(fnName{1});
      
      [testOut{1:length(expectedOut)}] = feval(fh, args.argin{:});
      
      % testOut fig handles
      % TODO: generalize fig handle checking
      if all(isValidFigHandle(testOut{1}))
        handles = testOut{1};
        
        dsUnitSave_figHandles( handles, tempDirPath )
        
        close all
      end
      
      % test output args
      for ind = 1:length(expectedOut)
        testCase.verifyEqual(testOut{ind}, expectedOut{ind});
      end
      
      % test output files
      expectedOutputDir = fullfile(dataDirPath, 'output');
      if exist(expectedOutputDir, 'dir')
        expectedOutputFiles = rls(expectedOutputDir);
        expectedOutputFiles(1) = []; % remove 'output' reference
        testOutputFiles = rls(tempDirPath);
        testOutputFiles(1) = []; % remove parent folder reference
        
        % test same num files
        testCase.assertLength(testOutputFiles, length(expectedOutputFiles));
        
        % test each file
        for iFile = 1:length(expectedOutputFiles)
%           thisExpectedOutputFilePath = fullfile(expectedOutputDir, expectedOutputFiles{iFile});
          thisExpectedOutputFilePath = expectedOutputFiles{iFile};
          
%           thistestOutputFilePath = fullfile(tempDirPath, testOutputFiles{iFile});
          thistestOutputFilePath = testOutputFiles{iFile};
          
          [~,~,ext] = fileparts2(thistestOutputFilePath);
          switch ext
            case '' % dir
              continue
            case '.fig'
              % compareFigFiles 
              testCase.verifyEmpty( compareFigFiles2(thistestOutputFilePath, thisExpectedOutputFilePath, true) );
            case '.mat'
              % compare loaded structs
              thistest = load(thistestOutputFilePath);
              thisExpected = load(thisExpectedOutputFilePath);
              testCase.verifyEqual(thistest, thisExpected);
            case '.mex4unittest'
              continue % skip testing of mex files
            otherwise
              % do diff
              testCase.verifyFalse( logical(system(sprintf('diff -qr %s %s', thistestOutputFilePath, thisExpectedOutputFilePath))) );
          end
        end
      end
      
    end
  end
  
end
