function saveAutoGenTestData(argin, argout)
% Inputs:
%   argin: cell array with input arguments
%   argout: cell array with output arguments

stack = dbstack;
fnName = stack(end-1).name;

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
