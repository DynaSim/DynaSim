function dsUnitSaveAutoGenTestDir(argin, argout, dirIn, dirOut, localFn_flag)
% Inputs:
%   argin: cell array with input arguments
%   argout: cell array with output arguments
%   dirin: path to directory that is input to fxn
%   localFn_flag: logical of whether to include 'fn1__fn2' in fnName or 'fn1'

if ~exist('localFn_flag','var')
  localFn_flag = false;
end

% get fn name from stack
stack = dbstack;
if ~localFn_flag
  fnNameStack = stack(2).name;
  fnName = fnNameStack;
else
  fnNameStack = stack(4).name;
  fnName = [fnNameStack '__' stack(3).name];
end

% check package namespace
if ~isfunction(fnNameStack) || strcmp(fnNameStack, 'strrep')
  if isfunction(['ds.' fnNameStack])
    fnName = ['ds.' fnName];
  elseif isfunction(['dsUnit.' fnNameStack])
    fnName = ['dsUnit.' fnName];
  end
end

hash = DataHash(argin);

% test dir
testDirName = sprintf('%s_autogen_%s', fnName, hash);
testFileDir = fullfile(dsGetConfig('ds_unitTestData_path'), 'autogenDirs_newSave', testDirName);
mkdirSilent(testFileDir);

% args file
testFileName = 'args.mat';
testFilePath = fullfile(testFileDir, testFileName);
save(testFilePath, 'argin', 'argout')

% output dir
if exist('dirOut', 'var') && ~isempty(dirOut)
  testOutputDir = fullfile(testFileDir, 'output');
  mkdirSilent(testOutputDir);
  
  studyDirInd = find(strcmp(argin, 'study_dir'));
  if ~isempty(studyDirInd) % if study_dir defined in argin
    studyDir = argin{studyDirInd+1};
    testOutputStudyDir = fullfile(testOutputDir,studyDir);
    mkdirSilent(testOutputStudyDir);
    
    copyfile(dirOut, testOutputStudyDir);
  else
    copyfile(dirOut, testOutputDir);
  end
end

% argout fig handles to output dir
% 	TODO: generalize fig handle checking
if all(isValidFigHandle(argout{1}))
  handles = argout{1};
  testOutputDir = fullfile(testFileDir, 'output');
  mkdirSilent(testOutputDir);
  
  dsUnitSave_figHandles( handles, testOutputDir )
end

% input dir
if exist('dirIn', 'var') && ~isempty(dirIn)
  testInputDir = fullfile(testFileDir, 'input');
  mkdirSilent(testInputDir);
  
  copyfile(dirIn, testInputDir);
end

end
