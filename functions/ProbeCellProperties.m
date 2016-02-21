function data = ProbeCellProperties(model,varargin)
% data = ProbeCellProperties(model,'option1',option1,...)
% This experiment is designed to probe the intrinsic properties of cells in
% one or more populations. It removes all connections between cells and 
% populations and then runs a series of simulations delivering
% hyperpolarizing and depolarizing pulses. It is designed to be used in
% conjunction with the analysis function "CalcCellProperties" which accepts
% the data array produced by this experiment and returns the
% electrophysiological properties characterizing each cell's response in
% the population.
% 
% Example: ...
% data=ProbeCellProperties(model,'verbose_flag',1);
% 
% Note: this function is based on the DNSim experiment "cell_pulses".
% See also: CalcCellProperties

% Experiment: input model, produces data sets for all step levels
% Analysis: input data sets for all step levels, output one stat structure  
%           per experiment call with ephys properties for each cell in each 
%           population of the model.

% Check inputs
options=CheckOptions(varargin,{...
  'target_equation','ODE1',[],...
  'amplitudes',0:2:10,[],...
  'tspan',[0 250],[],...
  'onset',50,[],...
  'offset',240,[],...
  },false);

model=CheckModel(model);

% Remove connections from the model specification and regenerate the model
if ~isempty(model.specification.connections)
  specification=model.specification;
  specification.connections=[];
  model=GenerateModel(specification);
end

% Extract population info
num_pops=length(model.specification.populations);
pop_names={model.specification.populations.name};

% Prepare list of modifications to add input pulses
modifications={};
for i=1:num_pops
  modifications(end+1,:)={pop_names{i},'equations',['cat(' options.target_equation ',+pulse(t); pulse(t)=TONIC*(t>=onset&t<=offset); monitor pulse)']};
  modifications(end+1,:)={pop_names{i},'onset',options.onset};
  modifications(end+1,:)={pop_names{i},'offset',options.offset};
end

% Prepare 'vary' specification to adjust pulse amplitudes in all populations
% simultaneously: {'(pop1,pop2,...)','TONIC',amplitudes}
objects='(';
for i=1:num_pops
  objects=[objects pop_names{i} ','];
end
objects=[objects(1:end-1) ')'];
vary={objects,'TONIC',options.amplitudes};  

% apply modifications to effectively add experimental apparatus to model
model=ApplyModifications(model,modifications);

% execute experimental protocol by varying parameters across simulations
fprintf('Running experiment: %s\n',mfilename);
keyvals=RemoveKeyval(varargin,'tspan');
data=SimulateModel(model,'vary',vary,'tspan',options.tspan,keyvals{:});

% add options to data
for i=1:length(data)
  data(i).model.experiment_options=options;
end
 
  