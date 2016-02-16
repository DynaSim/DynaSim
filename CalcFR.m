function data = CalcFR(data,varargin)
%% data = CalcFR(data,'option',value)
% Inputs:
%   data - DynaSim data structure (see CheckData)
%   options:
%     'variable' - name of field containing data on which to calculate firing
%                rates (default: *_spikes or first variable in data.labels)
%     'threshold' - scalar threshold value for detecting events (default: 0)
%     'bin_size' - size of temporal window over which to calculate rate [ms or fraction of data set] (default: 5% of the data set)
%     'bin_shift' - how much to shift the bin before calculating rate again [ms or fraction of data set] (default: 1% of the data set)
%     'exclude_data_flag' - whether to remove simulated data from result structure (default: 0)
% Outputs:
%   data: data structure with firing rates [Hz] in .variable_FR
% 
% note:
% "variable" can be specified as the name of a variable listed in
% data.labels, a cell array of string listing variable names, or as a 
% regular expression pattern for identifying variables to process.
% See SelectVariables for more info on supported specifications.
% 
% note:
% DynaSim spike monitor returns spike data in variables *_spikes.
% eg, data=SimulateModel('dv/dt=@current+10; {iNa,iK}; monitor v.spikes');
%     returns spikes in data.pop1_v_spikes (where 'pop1' is the default
%     population name if not specified by the user).
% 
% Examples:
% s.populations(1).name='E';
% s.populations(1).equations='dv/dt=@current+10; {iNa,iK}; v(0)=-65';
% s.populations(2).name='I';
% s.populations(2).equations='dv/dt=@current+10; {iNa,iK}; v(0)=-65';
% data=SimulateModel(s);
% data=CalcFR(data,'variable','*_v');
% data % contains firing rates for E and I pops in .E_v_FR and .I_v_FR.
% 
% See also: PlotFR, AnalyzeStudy, SimulateModel, CheckData, SelectVariables

%% 1.0 Check inputs
options=CheckOptions(varargin,{...
  'variable',[],[],...        
  'threshold',1e-5,[],... % slightly above zero in case variable is point process *_spikes {0,1}
  'bin_size',.05,[],...  % 30
  'bin_shift',.01,[],... % 10
  'exclude_data_flag',0,{0,1},...
  },false);

data = CheckData(data);
% note: calling CheckData() at beginning enables analysis function to
% accept data matrix [time x cells] in addition to DynaSim data structure.

if numel(data)>1
  % use AnalyzeStudy to recursively call CalcFR on each data set
  data=AnalyzeStudy(data,@CalcFR,varargin{:});
  return;
end

% time info
time = data.time;
dt = time(2)-time(1);
ntime=length(time);

% set defaults
% default variable to process
if isempty(options.variable)
  if any(~cellfun(@isempty,regexp(data.labels,'_spikes$')))
    % use results from DynaSim spike monitor
    options.variable=data.labels(~cellfun(@isempty,regexp(data.labels,'_spikes$')));
    if length(options.variable)==1 % store in string
      options.variable=options.variable{1};
    end
  else
    % use first state variable in model
    %options.variable=data.labels{1};
  end
end
% check bin_size
if options.bin_size>1
  % convert from ms to time points
  options.bin_size=ceil(options.bin_size/dt);
else
  % convert from fraction to time points
  options.bin_size=ceil(options.bin_size*ntime);
end
% constrain bin_size to entire data set
if options.bin_size>ntime
  options.bin_size=ntime;
end
% check bin_shift
if options.bin_shift>1
  % convert from ms to time points
  options.bin_shift=ceil(options.bin_shift/dt);
else
  % convert from fraction to time points
  options.bin_shift=ceil(options.bin_shift*ntime);
end

%% 2.0 set list of variables to process as cell array of strings
options.variable=SelectVariables(data(1).labels,options.variable);

%% 3.0 calculate firing rates for each variable
if ~isfield(data,'results')
  data.results={};
end
% 3.1 calc bin info
% samples at which bins begin
bin_index_begs=1:options.bin_shift:ntime;
% samples at which bins end
bin_index_ends=bin_index_begs+options.bin_size;
% remove final bin if extends beyond data
if bin_index_ends(end)>ntime
  bin_index_begs=bin_index_begs(bin_index_ends<=ntime);
  bin_index_ends=bin_index_ends(bin_index_ends<=ntime);
end
% times at which bins begin
bin_times=time(bin_index_begs);
% number of bins
nbins=length(bin_index_begs);
% time width of a single bin in seconds
bin_width=(dt/1000)*options.bin_size;
% 3.2 loop over variables to process
for v=1:length(options.variable)
  % extract this data set
  var=options.variable{v};
  dat=data.(var);
  % determine how many cells are in this data set
  ncells=size(dat,2);
  % loop over cells
  FR=zeros(nbins,ncells);
  spike_times=cell(1,ncells);
  for i=1:ncells
    % get spikes in this cell
    spike_inds=1+find((dat(2:end,i)>=options.threshold & dat(1:end-1,i)<options.threshold));
    spikes=zeros(ntime,1);
    spike_times{i}=time(spike_inds);
    if any(spike_inds)
      spikes(spike_inds)=1;
      % calculate firing rates
      for bin=1:nbins
        % (# spikes in bin) / (duration of bin in seconds)
        FR(bin,i)=sum(spikes(bin_index_begs(bin):bin_index_ends(bin)))/bin_width;
      end
    end
  end
  % add firing rates to data structure
  data.([var '_FR'])=FR;
  data.([var '_spike_times'])=spike_times;
  if ~ismember([var '_FR'],data.results)  
    data.results{end+1}=[var '_FR'];
    data.results{end+1}=[var '_spike_times'];
  end
end
% add bin times to data
data.time_FR=bin_times;
if ~ismember('time_FR',data.results)
  data.results{end+1}='time_FR';
end
if options.exclude_data_flag
  for l=1:length(data.labels)
    data=rmfield(data,data.labels{l});
  end
end

