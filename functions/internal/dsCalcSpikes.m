function data = dsCalcSpikes(data, varargin)
%dsCalcSpikes - Calculate spike indicies for DynaSim data
%
% Usage:
%   data = dsCalcSpikes(data,'option',value)
%
% Inputs:
%   - data: DynaSim data structure (see dsCheckData)
%   - options:
%     'variable'         : name of field containing data on which to calculate
%                          spike indicies (default: first variable in data.labels)
%     'threshold'        : scalar threshold value for detecting events (default: 0)
%     'overwrite_flag'   : whether to overwrite existing spikes fields (default: 0)
%     'output_suffix'    : suffix to add to result field name (default: '', i.e. none)
%     'exclude_data_flag': whether to remove simulated data from result
%                          structure (default: 0)
%
% Outputs:
%   - data: data structure with spike indicies in .variable_spikes
%
% Notes:
% - "variable" can be specified as the name of a variable listed in
%     data.labels, a cell array of string listing variable names, or as a regular
%     expression pattern for identifying variables to process. See dsSelectVariables
%     for more info on supported specifications.
% - DynaSim spike monitor returns spike data in variables *_spikes.
%   - e.g., `data=dsSimulate('dv/dt=@current+10; {iNa,iK}; monitor v.spikes');`
%     returns spikes in data.pop1_v_spikes (where 'pop1' is the default
%     population name if not specified by the user).
%
% Examples:
%   s.populations(1).name='E';
%   s.populations(1).equations='dv/dt=@current+10; {iNa,iK}; v(0)=-65';
%   s.populations(2).name='I';
%   s.populations(2).equations='dv/dt=@current+10; {iNa,iK}; v(0)=-65';
%   data = dsSimulate(s);
%   data = dsCalcSpikes(data,'variable','*_v');
%   data % contains spikes for E and I pops in .E_v_spikes and .I_v_spikes.
%
% See also: dsAnalyzeStudy, dsSimulate, dsCheckData, dsSelectVariables
% 
% Author: Erik Roberts
% Copyright (C) 2018 Jason Sherfey, Boston University, USA

%% 1.0 Check inputs
options=dsCheckOptions(varargin,{...
  'variable','',[],...
  'time_limits',[-inf inf],[],...
  'threshold',1e-5,[],... % slightly above zero in case variable is point process *_spikes {0,1}
  'overwrite_flag',0,{0,1},... % whether to overwrite existing spikes fields
  'exclude_data_flag',0,{0,1},...
  'output_suffix','',[],...
  'auto_gen_test_data_flag',0,{0,1},...
  },false);

%% auto_gen_test_data_flag argin
if options.auto_gen_test_data_flag
  varargs = varargin;
  varargs{find(strcmp(varargs, 'auto_gen_test_data_flag'))+1} = 0;
  varargs(end+1:end+2) = {'unit_test_flag',1};
  argin = [{data}, varargs]; % specific to this function
end

data = dsCheckData(data, varargin{:});
% note: calling dsCheckData() at beginning enables analysis function to
% accept data matrix [time x cells] in addition to DynaSim data structure.

if numel(data)>1
  % use dsAnalyzeStudy to recursively call dsCalcSpikes on each data set
  data = dsAnalyzeStudy(data,@dsCalcSpikes,varargin{:});
  return;
end

% time info
time = data.time;
% dt = time(2)-time(1);
t1 = nearest(time, options.time_limits(1)); % index to first sample
t2 = nearest(time, options.time_limits(2)); % index to last sample
% nTime = (t2-t1)+1;
timeInds = t1:t2;


%% 2.0 set list of variables to process as cell array of strings
options.variable = dsSelectVariables(data(1),options.variable, varargin{:});


%% 3.0 calculate spikes for each variable
if ~isfield(data,'results')
  data.results = {};
end

% 3.1 loop over variables to process
for iVar = 1:length(options.variable)
  % get variable name
  var = options.variable{iVar};
  thisVarStr = [var '_spikes' options.output_suffix];
  
  if isfield(data, thisVarStr) && ~options.overwrite_flag
    continue
  end
  
  % extract this data set
  dat = data.(var);
  
  % trim data
  dat = dat(timeInds, :,:);
  
  % determine how many cells are in this data set
  nTime = size(dat, 1);
  nCells = size(dat, 2);
  
  % loop over cells
  spikeProcess = zeros(nTime, nCells, 'single');
  for iCell = 1:nCells
    % get spikes in this cell ([0,1] increment process)
    % check for first threshold crossing
    spikeProcess(2:end,iCell) = (dat(2:end,iCell) >= options.threshold & dat(1:end-1,iCell) < options.threshold);
  end
  
  % add new spikes to data structure
  data.(thisVarStr) = spikeProcess;
  
  % check results field for var
  if ~ismember(thisVarStr, data.results)
    data.results{end+1} = thisVarStr;
  end
  
  % check labels field for var
  if ~ismember(thisVarStr, data.labels)
    data.labels{end+1} = thisVarStr;
  end
end

% check for removing data
if options.exclude_data_flag
  for l = 1:length(data.labels)
    thisLabel = data.labels{l};
    
    % only remove non spike vars
    if isempty(regexp(thisLabel,'_spikes$', 'once'))
      data = rmfield(data, thisLabel);
    end
  end
end

%% auto_gen_test_data_flag argout
if options.auto_gen_test_data_flag
  argout = {data}; % specific to this function
  
  dsUnitSaveAutoGenTestData(argin, argout);
end
