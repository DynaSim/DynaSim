function [ files ] = getAutogenFiles( localfn_flag )
%GETAUTOGENFILES returns cell list of autogen files
%
% Inputs:
%   localfn_flag: if true, only returns local function files. if false, does not
%                 return local fn files.

files = lscell(fullfile(ds.getConfig('ds_testData_path'), 'autogen', '*_autogen_*'));

if localfn_flag
  files = files{~cellfun(@isempty,regexp(files,'__'))};
else
  files = files{cellfun(@isempty,regexp(files,'__'))};
end

end

