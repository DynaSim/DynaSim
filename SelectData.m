function data=SelectData(data,varargin)
%% data=SelectData(data,'option',value)
% Purpose: select subset of data
% Inputs:
%   data -- DNSim data structure (see SimulateModel)
%   options:
%     'time_limits' -- [beg,end] (units of data.time)
% 
% Example:
% data=DataSelection(data,'time_limits',[20 80]); % return simulated data from 20-80ms
% 
% See also: SimulateModel, ImportData

% todo: specify subsets to return in terms of varied parameters, toilim, ROIs, etc
% possible format for specifying range_varied:
% {'E','gNa',[.1 .3]; 'I->E','tauI',[15 25]; 'I','mechanism_list','+iM'}

% check inputs
data=CheckData(data);

% recursively call SelectData if more than one data set
if length(data)>1
  for i=1:length(data)
    data(i)=SelectData(data(i),varargin{:});
  end
  return;
end

options=CheckOptions(varargin,{...
  'time_limits',[-inf inf],[],...
  },false);

% select subset from state variables and monitors
time=data.time;
seltime=time>=options.time_limits(1) & time<=options.time_limits(2);
for i=1:length(data.labels)
  data.(data.labels{i})=data.(data.labels{i})(seltime,:);
end
  
% todo: select subset from time-varying post-processed results
if isfield(data,'results')
  % get time vectors: time_(*) with * matching matching names in data.results
  % ...
  
end
  
