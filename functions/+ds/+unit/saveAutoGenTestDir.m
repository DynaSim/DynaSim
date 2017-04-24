function saveAutoGenTestDir(argin, argout, dirin, localFn_flag)
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
if ~isfunction(fnNameStack)
  if isfunction(['ds.' fnNameStack])
    fnName = ['ds.' fnName];
  end
end

hash = DataHash(argin);

% test dir
testDirName = sprintf('%s_autogen_%s', fnName, hash);
testFileDir = fullfile(ds.getConfig('ds_testData_path'), 'autogenDirs_newSave', testDirName);
mkdirSilent(testFileDir);

% argout fig handles
% TODO: generalize fig handle checking
if all(isValidFigHandle(argout{1}))
  handles = argout{1};
  testOutputDir = fullfile(testFileDir, 'output');
  mkdirSilent(testOutputDir);
  
  ds.unit.save_figHandles( handles, testOutputDir )
end

% args file
testFileName = 'args.mat';
testFilePath = fullfile(testFileDir, testFileName);
save(testFilePath, 'argin', 'argout')

% input dir
if exist('dirin', 'var')
  testInputDir = fullfile(testFileDir, 'input');
  mkdirSilent(testInputDir);
  copyfile(dirin, testInputDir)
end

end
