function saveAutoGenTestData(argin, argout, localFn_flag)
% Inputs:
%   argin: cell array with input arguments
%   argout: cell array with output arguments
%   localFn_flag: logical of whether to include 'fn1__fn2' in fnName or 'fn1'

if ~exist('localFn_flag','var')
  localFn_flag = false;
end

stack = dbstack;
fnName = stack(end-1).name;

if localFn_flag
  fnName = [fnName '__' stack(end-2).name]
end

% check package namespace
if ~isfunction(fnName)
  if isfunction(['ds.' fnName])
    fnName = ['ds.' fnName];
  end
end

hash = DataHash(argin);

testFileName = sprintf('%s_autogen_%s.mat', fnName, hash);
testFileDir = fullfile(ds.getConfig('ds_testData_path'), 'autogen');
mkdirSilent(testFileDir);
testFilePath = fullfile(testFileDir, testFileName);
save(testFilePath, 'argin', 'argout')

end
