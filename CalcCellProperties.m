function stats = CalcCellProperties(data,varargin)
% stats = CalcCellProperties(data,'option1',option1,...)
% This analysis function calculates the intrinsic electrophysiological 
% properties of all cells in one or more populations. It is designed to be 
% used in conjunction with the experiment "ProbeCellProperties" which 
% removes all connections from a model and produces a data array of 
% simulated data in response to a series of hyperpolarizing and 
% depolarizing pulses.
% 
% Inputs:
% - data -- array of DynaSim data structures returned by ProbeCellProperties
% 
% Properties measured:
% - RMP: resting membrane potential [mV]
%       Method: ... (see [1])
% - ...
% 
% References for methods used:
% [1] ...
% 
% Example: ...
% data = ProbeCellProperties(model)
% stats = CalcCellProperties(data)
%
% Note: this function is based on the DNSim experiment "cell_pulses".
% See also: ProbeCellProperties

% Check inputs
options=CheckOptions(varargin,{...
  'spike_detection_threshold',.5,[],...
  'skip_time',10,[],... % time [ms] to exclude from detection
  'plot_flag',1,{0,1},...
  },false);

data=CheckData(data);
model=data(1).model;
time=data(1).time;
num_times=length(time);

% extract population info
pop_names={model.specification.populations.name};
pop_sizes=[model.specification.populations.size];
num_pops=length(pop_names);

% get list of variables to analyze
[vars,vars_pop]=SelectVariables(data(1).labels);

% assume only one variable per population
if length(unique(vars_pop))>num_pops
  [vars;vars_pop]
  error('there can only be one variable per population.');
end

% extract experiment parameters
experiment_options=data(1).model.experiment_options;

% stimulus amplitude
amplitudes=experiment_options.amplitudes; %[data.(data(1).varied{1})];
num_amps=length(amplitudes);

% stimulus timing
onset=experiment_options.onset; %model.parameters.([pop_names{1} '_onset']);
offset=experiment_options.offset; %model.parameters.([pop_names{1} '_offset']);
tsel=time>=onset&time<=offset;

% assume exactly one simulation per amplitude
if num_amps ~= length(data)
  error('there can only be one simulation per stimulus amplitude');
end

% initialize stats structure
% analyze each population
stats.amplitudes=amplitudes;
stats.time=time;
stats.cell_results={}; % list of field names storing results for each cell
stats.pop_results={};  % list of field names storing results for each population

% collect pulses
input=zeros(num_times,num_amps);
for a=1:num_amps
  input(:,a)=data(a).([pop_names{1} '_pulse'])(:,1);
end

% analyze each population
for p=1:num_pops
  % extract info for this population
  this_pop=pop_names{p};
  this_var=vars{strcmp(vars_pop,this_pop)};
  num_cells=pop_sizes(p);
  % collect data for this population across simulations
  % note: each simulation has a different input amplitude
  X=zeros(num_amps,num_times,num_cells);
  for a=1:num_amps
    X(a,:,:)=data(a).(this_var);
  end
  
  % calculate simple means
  means=nanmean(X(:,tsel,:),2); % average over select times
  % store means in stats structure
  stats.([this_var '_mean_per_amp'])=squeeze(means);
  stats.([this_var '_pop_mean_per_amp'])=nanmean(means,3);
  stats.([this_var '_pop_mean_per_time'])=nanmean(X,3);
%   stats.cell_results{end+1}=[this_var '_mean'];
%   rstats.pop_results{end+1}=[this_var '_pop_mean'];
  
  % calculate more sophisticated measures for each cell
  for c=1:num_cells
    % ...
  end
end

if options.plot_flag
  v=1;
  figure('position',[180 330 1450 530])
  subplot(2,2,1); % V(t)
  plot(time,stats.([vars{v} '_pop_mean_per_time']));
  xlabel('time [ms]'); ylabel(['<' strrep(vars{v},'_','\_') ', pop>']);
  legend(cellfun(@(x)['I=' num2str(x)],num2cell(amplitudes),'uni',0));
  subplot(2,2,3); % I(t)
  plot(time,input); xlabel('time [ms]'); ylabel('I(t)');
  legend(cellfun(@(x)['I=' num2str(x)],num2cell(amplitudes),'uni',0));
  subplot(2,2,2); % I/V
  plot(amplitudes,stats.([vars{v} '_mean_per_amp'])); hold on
  plot(amplitudes,stats.([vars{v} '_pop_mean_per_amp']),'k-','linewidth',5);
  xlabel('amplitude'); ylabel(['<' strrep(vars{v},'_','\_') ', time>']);
end


  
  
  
  
  
  