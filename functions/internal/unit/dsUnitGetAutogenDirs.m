function [ files ] = dsUnitGetAutogenDirs( localfn_flag, query_flag )
%GETAUTOGENFILES returns cell list of autogen directories, for functions which write files to disk
%
% Inputs:
%   localfn_flag: if true, only returns local function files. if false, does not
%                 return local fn files. Default: false;
%   query_flag: whether to use query variable from base

if ~exist('localfn_flag','var')
  localfn_flag = false;
end
if ~exist('query_flag','var')
  query_flag = false;
end

if ~query_flag
  files = lscell(fullfile(dsGetConfig('ds_unitTestData_path'), 'autogenDirs'));
else
  % TODO: fix this
  files = lscell(fullfile(dsGetConfig('ds_unitTestData_path'), 'autogenDirs', ['*' evalin('base','query') '*_autogen_*']));
end

if localfn_flag
  files = files( ~cellfun(@isempty,regexp(files,'__')) );
else
  files = files( cellfun(@isempty,regexp(files,'__')) );
end

end
