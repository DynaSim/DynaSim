function [param_mat,varied,param_cell]=CollectVariedParams(data)
% [all_values,param_names,unique_values]=CollectVariedParams(data)
% inputs:
%   data: 
% outputs:
%   all_values: [num_sets x num_params_varied], values used for each data set
%   param_names: list of names of varied parameters
%   unique_values: cell array of unique values used for each varied  parameter

% Check inputs
if ~isstruct(data)
  fprintf('input must be a structure...exiting CollectVariedParams.\n');
  return;
end
if ~isfield(data,'varied')
  fprintf('no varied info in data.\n');
  return;
end

% collect info on parameters varied
varied=data(1).varied;
num_varied=length(varied); % number of model components varied across simulations
num_sims=length(data); % number of data sets (one per simulation)
% collect info on parameters varied
param_mat=zeros(num_sims,num_varied); % values for each simulation
param_cell=cell(1,num_varied); % unique values for each parameter
% loop over varied components and collect values
for j=1:num_varied
  if isnumeric(data(1).(varied{j}))
    param_mat(:,j)=[data.(varied{j})]; % values for each simulation
    param_cell{j}=unique([data.(varied{j})]); % unique values for each parameter
  else
    % todo: handle sims varying non-numeric model components 
    % (eg, mechanisms) (also in PlotFR and SelectData)
  end
end
