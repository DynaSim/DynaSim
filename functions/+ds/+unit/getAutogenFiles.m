function [ files ] = getAutogenFiles( localfn_flag, query_flag )
%GETAUTOGENFILES returns cell list of autogen files
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
  files = lscell(fullfile(ds.getConfig('ds_testData_path'), 'autogen', '*_autogen_*'));
else
  files = lscell(fullfile(ds.getConfig('ds_testData_path'), 'autogen', ['*' evalin('base','query') '*_autogen_*']));
end

if localfn_flag
  files = files( ~cellfun(@isempty,regexp(files,'__')) );
else
  files = files( cellfun(@isempty,regexp(files,'__')) );
end

end

