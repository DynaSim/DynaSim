function dsUnitSaveAutoGenTestData(argin, argout, localFn_flag)
% Inputs:
%   argin: cell array with input arguments
%   argout: cell array with output arguments
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

testFileName = sprintf('%s_autogen_%s.mat', fnName, hash);
testFileDir = fullfile(dsGetConfig('ds_unitTestData_path'), 'autogen_newSave');
mkdirSilent(testFileDir);
testFilePath = fullfile(testFileDir, testFileName);
save(testFilePath, 'argin', 'argout')

end
