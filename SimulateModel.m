function [data,studyinfo]=SimulateModel(model,varargin)
%% data=SimulateModel(model,'option',value,...)
% Purpose: manage simulation of a DynaSim model. This high-level function 
% offers many options to control the number of simulations, how the model is 
% optionally varied across simulations, and how the numerical integration
% is performed. It can optionally save the simulated data and create/submit
% simulation jobs to a compute cluster.
% Inputs:
%   model: DynaSim model structure or equations (see GenerateModel and 
%          CheckModel for more details)
% 
%   solver options (provided as key/value pairs: 'option1',value1,'option2',value2,...):
%     'solver'      : solver for numerical integration (see GetSolveFile)
%                     {'euler','rk2','rk4'} (default: 'rk4')
%     'tspan'       : time limits of simulation [begin,end] (default: [0 100]) [ms]
%                     note: units must be consistent with dt and model equations
%     'dt'          : time step used for DynaSim solvers (default: .01) [ms]
%     'downsample_factor': downsampling applied during simulation (default: 1, no downsampling) 
%                     (only every downsample_factor-time point is stored in memory and/or written to disk)
%     'ic'          : numeric array of initial conditions, one value per state 
%                     variable (default: all zeros). overrides definition in model structure
%     'random_seed' : seed for random number generator (default: 'shuffle', set randomly) (usage: rng(options.random_seed))
%     'compile_flag': whether to compile simulation using coder instead of 
%                     interpreting Matlab {0 or 1} (default: 0)
% 
%   options for running sets of simulations:
%     'vary'        : (default: [], vary nothing): cell matrix specifying model
%                     components to vary across simulations (see NOTE 1 and Vary2Modifications)
% 
%   options to control saved data:
%     'save_results_flag': whether to save results of analysis and plotting
%     'save_data_flag': whether to save simulated data to disk after completion {0 or 1} (default: 0)
%     'overwrite_flag': whether to overwrite existing data files {0 or 1} (default: 0)
%     'study_dir'     : relative or absolute path to output directory (default: current directory)
%     'prefix'        : string to prepend to all output file names (default: 'study')
%     'disk_flag'     : whether to write to disk during simulation instead of storing in memory {0 or 1} (default: 0)
%     'precision'     : {'single','double'} precision of simulated data saved to disk (default: 'single')
%                 
%   options for cluster computing:
%     'cluster_flag'  : whether to run simulations on a cluster submitted 
%                     using qsub (see CreateBatch) {0 or 1} (default: 0)
%     'sims_per_job'  : number of simulations to run per batch job (default: 1)
%     'memory_limit'  : memory to allocate per batch job (default: '8G')
% 
%   options for parallel computing: (requires Parallel Computing Toolbox)
%     'parallel_flag' : whether to use parfor to run simulations {0 or 1} (default: 0)
%     'num_cores'     : number of cores to specify in the parallel pool
%     *note: parallel computing has been disabled for debugging...
% 
%   options for post-processing:
%     'analysis_functions': cell array of analysis function handles
%     'analysis_options'  : cell array of option cell arrays {'option1',value1,...}
%     'plot_functions'    : cell array of plot function handles
%     'plot_options'      : cell array of option cell arrays {'option1',value1,...}
% 
%   other options:
%     'verbose_flag'  : whether to display informative messages/logs (default: 0)
%     'modifications' : how to modify DynaSim specification structure component before simulation (see ApplyModifications)
%     'experiment'    : function handle of experiment function (see NOTE 2)
%     'experiment_options' : single cell array of key/value options for experiment function
%     'optimization'  : function handle of optimization function (see NOTE 2)
% 
% Outputs:
%   DynaSim data structure:
%     data.labels           : list of state variables and monitors recorded
%     data.(state_variables): state variable data matrix [time x cells]
%     data.(monitors)       : monitor data matrix [time x cells]
%     data.time             : time vector [time x 1]
%     data.simulator_options: simulator options used to generate simulated data
%     data.model            : model used to generate simulated data
%     [data.varied]         : list of varied model components (present only if anything was varied)
% 
%   DynaSim studyinfo structure (only showing select fields, see CheckStudyinfo for more details)
%     studyinfo.study_dir
%     studyinfo.base_model (=[]): original model from which a set of simulations was derived
%     studyinfo.base_simulator_options (=[])
%     studyinfo.base_solve_file (='')
%     studyinfo.simulations(k): metadata for each simulation in a set of simulations
%                           .sim_id         : unique identifier in study
%                           .modifications  : modifications made to the base model during this simulation
%                           .data_file      : full filename of eventual output file
%                           .batch_dir (=[]): directory where batch jobs were saved (if cluster_flag=1)
%                           .job_file (=[]) : m-file batch job that runs this simulation (if cluster_flag=1)
%                           .simulator_options: simulator options for this simulation
%                           .solve_file     : full filename of m- or mex-file that numerically integrated the model
% 
% NOTE 1: 'vary' indicates the variable to vary, the values
% it should take, and the object whose variable should be varied. 
% Syntax: vary={object, variable, values; ...}. For instance, to vary
% parameter 'gNa', taking on values 100 and 120, in population 'E', set
% vary={'E','gNa',[100 120]}. To additionally vary 'gSYN' in the connection
% mechanism from 'E' to 'I', set vary={'E','gNa',[100 120];'E->I','gSYN',[0 1]}.
% Mechanism lists and equations can also be varied. (see Vary2Modifications 
% for more details and examples).
% 
% EXAMPLES:
% Example 1: Lorenz equations with phase plot
%   eqns={
%     's=10; r=27; b=2.666';
%     'dx/dt=s*(y-x)';
%     'dy/dt=r*x-y-x*z';
%     'dz/dt=-b*z+x*y';
%   };
%   data=SimulateModel(eqns,'tspan',[0 100],'ic',[1 2 .5]);
%   plot(data.pop1_x,data.pop1_z); title('Lorenz equations'); xlabel('x'); ylabel('z')
% 
% Example 2: Leaky integrate-and-fire with spike monitor
%   eqns={
%     'tau=10; R=10; E=-70; I=1.55; thresh=-55; reset=-75';
%     'dV/dt=(E-V+R*I)/tau; if(V>thresh)(V=reset)';
%     'monitor V.spikes(thresh)';
%   };
%   data=SimulateModel(eqns,'tspan',[0 200],'ic',-75);
%   data.pop1_V(data.pop1_V_spikes==1)=20; % insert spike
%   plot(data.time,data.pop1_V); xlabel('time (ms)'); ylabel('V'); title('LIF with spikes')
% 
% Example 3: Hodgkin-Huxley-type Intrinsically Bursting neuron
%   eqns='dv/dt=5+@current; {iNaF,iKDR,iM}; gNaF=100; gKDR=5; gM=1.5; v(0)=-70';
%   data=SimulateModel(eqns,'tspan',[0 200]);
%   figure; plot(data.time,data.(data.labels{1}))
%   xlabel('time (ms)'); ylabel('membrane potential (mV)'); title('Intrinsically Bursting neuron')
% 
% Example 4: varying max Na+ conductance in Hodgkin-Huxley neuron
%   eqns='dv/dt=@current+10; {iNa,iK}; v(0)=-60';
%   data=SimulateModel(eqns,'vary',{'','gNa',[50 100 200]});
%   % plot how mean firing rate varies with parameter
%   PlotFR(data,'bin_size',30,'bin_shift',10); % bin_size and bin_shift in [ms]
% 
% Example 5: Sparse Pyramidal-Interneuron-Network-Gamma rhythm with rastergram
%   % define equations of cell model (same for E and I populations)
%   eqns={ 
%     'dv/dt=Iapp+@current/Cm+noise*randn(1,N_pop)*sqrt(dt)/dt';
%     'monitor v.spikes, iGABAa.functions, iAMPA.functions'
%   };
%   % define specification for two-population network model
%   s=[];
%   s.populations(1).name='E';
%   s.populations(1).size=80;
%   s.populations(1).equations=eqns;
%   s.populations(1).mechanism_list={'iNa','iK'};
%   s.populations(1).parameters={'Iapp',5,'gNa',120,'gK',36,'Cm',1,'noise',4};
%   s.populations(2).name='I';
%   s.populations(2).size=20;
%   s.populations(2).equations=eqns;
%   s.populations(2).mechanism_list={'iNa','iK'};
%   s.populations(2).parameters={'Iapp',0,'gNa',120,'gK',36,'Cm',1,'noise',4};
%   s.connections(1).source='I';
%   s.connections(1).target='E';
%   s.connections(1).mechanism_list={'iGABAa'};
%   s.connections(1).parameters={'tauD',10,'gSYN',.1,'netcon','ones(N_pre,N_post)'};
%   s.connections(2).source='E';
%   s.connections(2).target='I';
%   s.connections(2).mechanism_list={'iAMPA'};
%   s.connections(2).parameters={'tauD',2,'gSYN',.1,'netcon',ones(80,20)};
%   % simulate model
%   data=SimulateModel(s);
%   % plot voltages and rastergram
%   figure;
%   subplot(2,1,1); % voltage traces
%   plot(data.time,data.E_v,'b-',data.time,data.I_v,'r-')
%   title('Sparse Pyramidal-Interneuron-Network-Gamma (sPING)'); ylabel('membrane potential (mV)');
%   subplot(2,1,2); % rastergram
%   E_spikes=nan(size(data.E_v_spikes)); E_spikes(data.E_v_spikes==1)=1;
%   I_spikes=nan(size(data.I_v_spikes)); I_spikes(data.I_v_spikes==1)=1;
%   plot(data.time,E_spikes+repmat(1:80,[length(data.time) 1]),'bo'); hold on
%   plot(data.time,I_spikes+repmat(80+(1:20),[length(data.time) 1]),'ro'); axis([0 100 0 100]);
%   title('rastergram'); xlabel('time (ms)'); ylabel('cell index');
%   % simulate model varying two parameters (Iapp and tauD in sPING)
%   % warning: this may take up to a minute to complete:
%   vary={
%     'E'   ,'Iapp',[0 10 20];     % amplitude of tonic input to E-cells
%     'I->E','tauD',[5 10 15]      % inhibition decay time constant from I to E
%     };
%   data=SimulateModel(s,'vary',vary);
%   % plot firing rates calculated from spike monitor in both populations
%   PlotFR(data,'variable','*_spikes','bin_size',30,'bin_shift',10);
% 
% 
% See also: GenerateModel, CheckModel, GetSolveFile, CheckData, 
%           Vary2Modifications, CheckStudyinfo, CreateBatch

% todo: rename 'disk_flag' to something more descriptive

% dependencies: WriteDynaSimSolver, WriteMatlabSolver, PropagateFunctions, CheckModel,
% CheckOptions, Options2Keyval, DisplayError, DynaSim2Odefun

% <-- temporarily removed from help section -->
% NOTE 2: special functions that recursively call SimulateModel:
% - "Experiments" are ways of hacking the ODE system to incorporate additional 
% models (e.g., controlled inputs) and use them to simulate experimental 
% protocols by systematically varying Model Components across simulations in 
% prescribed ways. Technically, an Experiment could be any function that takes 
% a DynaSim model structure as its first input followed by key/value options. 
% Ideally, Experiments represent standardized procedural methods (experimental 
% protocols) for studying the modeled system. Experiment functions typically 
% involve applying a set of Modifications to a Base Model and varying the 
% modified model in prescribed ways. 
% - during "Optimization", each iteration involves a single Study producing 
% modified models, their simulated data sets and analysis results (e.g., 
% cost functions) that shape the Base Model for a subsequent iteration and 
% its Study. Hence, a closed-loop optimization protocol produces a set of 
% evolving Studies. Technically, an Optimization could be any function
% that takes a DynaSim model structure as its first input followed by
% key/value options. Optimization functions will typically involve
% while-looping through multiple Studies and analyzing data sets to update
% Model Components on each iteration until some stop condition is reached.

% Initialize outputs
data=[];
studyinfo=[];

% Check inputs
varargin = backward_compatibility(varargin);
options=CheckOptions(varargin,{...
  'tspan',[0 100],[],...          % [beg,end] (units must be consistent with dt and equations)
  'ic',[],[],...                  % initial conditions (overrides definition in model structure; can input as IC structure or numeric array)
  'solver','rk4',{'euler','rk1','rk2','rk4','modified_euler','rungekutta','rk','ode23','ode45'},... % DynaSim and built-in Matlab solvers
  'matlab_solver_options',[],[],... % options from odeset for use with built-in Matlab solvers
  'dt',.01,[],...                 % time step used for fixed step DynaSim solvers
  'downsample_factor',1,[],...    % downsampling applied during simulation (only every downsample_factor-time point is stored in memory or written to disk)
  'reduce_function_calls_flag',1,{0,1},...   % whether to eliminate internal (anonymous) function calls
  'save_parameters_flag',1,{0,1},...
  'random_seed','shuffle',[],...        % seed for random number generator (usage: rng(random_seed))
  'data_file','data.csv',[],... % name of data file if disk_flag=1
  'precision','single',{'single','double'},...
  'logfid',1,[],...
  'store_model_flag',1,{0,1},...  % whether to store model structure with data
  'verbose_flag',0,{0,1},...
  'modifications',[],[],...       % *DynaSim modifications structure
  'vary',[],[],...                % specification of things to vary or custom modifications_set
  'experiment',[],[],...          % experiment function. func(model,args)
  'experiment_options',[],[],...
  'optimization',[],[],...
  'cluster_flag',0,{0,1},...      % whether to run simulations on a cluster
  'sims_per_job',1,[],... % how many sims to run per batch job
  'memory_limit','8G',[],... % how much memory to allocate per batch job
  'parallel_flag',0,{0,1},...     % whether to run simulations in parallel (using parfor)
  'num_cores',4,[],... % # cores for parallel processing (SCC supports 1-12)
  'compile_flag',0,{0,1},... % exist('codegen')==6, whether to compile using coder instead of interpreting Matlab
  'disk_flag',0,{0,1},...            % whether to write to disk during simulation instead of storing in memory
  'save_data_flag',0,{0,1},...  % whether to save simulated data
  'save_results_flag',0,{0,1},...  % whether to save simulated data
  'project_dir',pwd,[],...
  'study_dir',[],[],... % study directory
  'prefix','study',[],... % prefix prepended to all output files
  'overwrite_flag',0,{0,1},... % whether to overwrite existing data
  'solve_file',[],[],... % m- or mex-file solving the system
  'sim_id',[],[],... % sim id in an existing study
  'studyinfo',[],[],... 
  'email',[],[],... % email to send notification upon study completion
  'analysis_functions',[],[],...
  'analysis_options',[],[],...
  'plot_functions',[],[],...
  'plot_options',[],[],...
  },false);
% more options: remove_solve_dir, remove_batch_dir, post_downsample_factor

if options.parallel_flag
  %error('parallel computing has been disabled for debugging. ''set parallel_flag'' to 0');
end

if options.compile_flag && options.reduce_function_calls_flag==0
  fprintf('setting ''reduce_function_calls_flag'' to 1 for compatibility with ''compile_flag''=1 (coder does not support anonymous functions).\n');
  options.reduce_function_calls_flag=1;
end

if options.cluster_flag && options.save_data_flag==0
%   options.save_data_flag=1;
%   if options.verbose_flag
%     fprintf('setting ''save_data_flag'' to 1 for storing results of batch jobs for later access.\n');
%   end
  options.save_results_flag=1;
  if options.verbose_flag
    fprintf('setting ''save_results_flag'' to 1 for storing results of batch jobs for later access.\n');
  end
end
% 
% if ischar(options.study_dir) && options.save_data_flag==0
%   options.save_data_flag=1;
%   if options.verbose_flag
%     fprintf('setting ''save_data_flag'' to 1 for storing results in study_dir: %s.\n',options.study_dir);
%   end
% end

% prepare analysis functions and options
if ~isempty(options.analysis_functions)
  if ~iscell(options.analysis_functions)
    % convert function handle into cell array of function handles
    options.analysis_functions={options.analysis_functions};
  end
  if any(~cellfun(@(x)isa(x,'function_handle'),options.analysis_functions))
    error('at least one analysis function was not provided as a function handle.');
  end
  if isempty(options.analysis_options)
    % convert to empty option cell array
    options.analysis_options={};
  end
  if ~iscell(options.analysis_options)
    error('''analysis_options'' must be a cell array of options or option cell arrays');
  end
  % force to be a cell array of option cell arrays
  if isempty(options.analysis_options) || ischar(options.analysis_options{1}) % first element is an option
    options.analysis_options={options.analysis_options};
  end
  % make sure there is one option cell array per analysis function
  if length(options.analysis_options)==1 && length(options.analysis_functions)>1
    % copy options for each analysis function
    options.analysis_options=repmat(options.analysis_options,[1 length(options.analysis_functions)]);
  elseif length(options.analysis_options) ~= length(options.analysis_functions)
    error('there must be one option cell array per analysis function.');
  end
%   if options.cluster_flag~=1
%     warning('analysis functions will not be run after simulation. currently automatic post-simulation analyses are supported only for cluster jobs.');
%     options.analysis_functions=[];
%     options.analysis_options=[];
%   end
end

% prepare plot functions and options
if ~isempty(options.plot_functions)
  if ~iscell(options.plot_functions)
    % convert function handle into cell array of function handles
    options.plot_functions={options.plot_functions};
  end
  if any(~cellfun(@(x)isa(x,'function_handle'),options.plot_functions))
    error('at least one plot function was not provided as a function handle.');
  end
  if isempty(options.plot_options)
    % convert to empty option cell array
    options.plot_options={};
  end
  if ~iscell(options.plot_options)
    error('''plot_options'' must be a cell array of options or option cell arrays');
  end
  % force to be a cell array of option cell arrays
  if isempty(options.plot_options) || ischar(options.plot_options{1}) % first element is an option
    options.plot_options={options.plot_options};
  end
  % make sure there is one option cell array per plot function
  if length(options.plot_options)==1 && length(options.plot_functions)>1
    % copy options for each plot function
    options.plot_options=repmat(options.plot_options,[1 length(options.plot_functions)]);
  elseif length(options.plot_options) ~= length(options.plot_functions)
    error('there must be one option cell array per plot function.');
  end
%   if options.cluster_flag~=1
%     warning('plot functions will not be run after simulation. currently automatic post-simulation plotting are supported only for cluster jobs.');
%     options.plot_functions=[];
%     options.plot_options=[];
%   end  
end

% check path
dynasim_path=fileparts(which(mfilename));
onPath=~isempty(strfind(path,[dynasim_path, pathsep]));
if ~onPath
  if options.verbose_flag
    fprintf('adding dynasim directory to Matlab path: %s\n',dynasim_path);
  end
  addpath(dynasim_path); % necessary b/c of changing directory for simulation
end
dynasim_functions=fullfile(dynasim_path,'functions');
onPath=~isempty(strfind(path,[dynasim_functions, pathsep]));
if ~onPath
  if options.verbose_flag
    fprintf('adding dynasim functions directory to Matlab path: %s\n',dynasim_functions);
  end
  addpath(dynasim_functions); % necessary b/c of changing directory for simulation
end

%% 1.0 prepare model and study structures for simulation

% handle special case of input equations with vary() statement
[model,options]=extract_vary_statement(model,options);

% check/standardize model
model=CheckModel(model); % handles conversion when input is a string w/ equations or a DynaSim specification structure

% 1.1 apply modifications before simulation and optional further variation across simulations
if ~isempty(options.modifications)
  [model,options.modifications]=ApplyModifications(model,options.modifications);
end

% 1.2 incorporate user-supplied initial conditions
if ~isempty(options.ic)
  if isstruct(options.ic) 
  % todo: create subfunc that converts numeric ICs to strings for
  % consistency (call here and in ProcessNumericICs <-- use code already there)
    % user provided structure with ICs
    warning('off','catstruct:DuplicatesFound');
    model.ICs=catstruct(model.ICs,options.ic);
  elseif isnumeric(options.ic)
    % user provided numeric array with one value per state variable per
    % cell with state variables ordered according to model.state_variables.    
    model.ICs=ProcessNumericICs;
  end
end

% expand set of things to vary across simulations
if ~isempty(options.vary)
  modifications_set=Vary2Modifications(options.vary,model);
else
  modifications_set={[]};
end

% 1.3 check for parallel simulations
% 1.3.1 manage cluster computing
% whether to write jobs for distributed processing on cluster
if options.cluster_flag==1
  % add to model any parameters in 'vary' not explicit in current model
  % approach: use ApplyModifications(), it does that automatically
  for i=1:length(modifications_set)
    if ~isempty(modifications_set{i}) && ~strcmp(modifications_set{i}{2},'mechanism_list') && ~strcmp(modifications_set{i}{2},'equations')
      model=ApplyModifications(model,modifications_set{i});
      break
    end
  end
  keyvals = Options2Keyval(options);
  studyinfo=CreateBatch(model,modifications_set,'simulator_options',options,'process_id',options.sim_id,keyvals{:});
  %if options.overwrite_flag==0
    % check status of study
%     [~,s]=MonitorStudy(studyinfo.study_dir,'verbose_flag',0,'process_id',options.sim_id);
%     if s==1 % study finished
%       if options.verbose_flag
%         fprintf('Study already finished. Importing data...\n');
%       end
%       studyinfo=CheckStudyinfo(studyinfo.study_dir,'process_id',options.sim_id);
%       data=ImportData(studyinfo,'process_id',options.sim_id);
%     end  
  %end
  return;
end

% 1.3.2 manage parallel computing on local machine 
    % TODO: debug local parallel sims, doesn't seem to be working right... 
    % (however SCC cluster+parallel works)
if options.parallel_flag==1
  % prepare solve_file
  tmp_options=options;
  if isempty(options.study_dir)
    tmp_options.study_dir=pwd;
  end
  if isempty(options.solve_file) || ~exist(options.solve_file,'file')
    solve_file=GetSolveFile(model,[],tmp_options);
  else
    solve_file=options.solve_file;
  end
  % prepare options
  keyvals=Options2Keyval(rmfield(options,{'vary','modifications','solve_file','parallel_flag'}));
  % open pool for distributed processing
  % or instead: require user to open pool before calling SimulateModel...
%   parpool(options.num_cores) 
  % run embarrassingly-parallel simulations
  
  
  % Previous parallel code would overwrite the same params.mat file on each
  % parallel iteration, resulting in the same parameters being used for all
  % simulations. This code circumvents this issue by assigning a different
  % study directory to each simulation.
  
  % List any core files - these should be deleted, as they are huge (debug)
%   system (['ls ' fullfile(options.study_dir,'output*')],'-echo');
%   system('find * -name "core*"','-echo');
 
  % Generate unique* ID based on date and time (*as long as they dont start on the
  % exact same second...). Using these unique identfiers will now enable
  % you to run multiple sims in the same folder.
  uniqueID = datestr(now,30);
  
  
  for sim = 1:length(modifications_set)
      mystudydirs{sim} = fullfile(options.study_dir,['output_parfor_' uniqueID '_' num2str(sim)]);                           
      
      % Create solve folders as needed
      if ~exist(mystudydirs{sim},'dir')                
          mkdir(fullfile(mystudydirs{sim}));           
          if ~exist(fullfile(mystudydirs{sim},'solve'),'dir') 
              mkdir(fullfile(mystudydirs{sim},'solve'));           
          end
      end
      
      
      [success,msg]=copyfile([strrep(solve_file,'_mex','') '*'],fullfile(mystudydirs{sim},'solve'));    % Copy the mex file into each study directory, to avoid re-compiling
      if ~success, error(msg); end
      
  end
 
  clear data
  parfor sim=1:length(modifications_set)
    %data(sim)=SimulateModel(model,'modifications',modifications_set{sim},'solve_file',solve_file,keyvals{:});       % Original parfor code
    data(sim)=SimulateModel(model,'modifications',modifications_set{sim},keyvals{:},'study_dir',mystudydirs{sim});  % My modification; now specifies a separate study directory for each sim
    disp(sim);
  end

  
% Clean up files leftover from sim
% Unfortunately we can't remove the folders due to locked .nfs files.
% Need to do this manually later...
  for sim = 1:length(mystudydirs)
      
    % Remove any "core" files that might be present and taking up space
    if exist(fullfile(mystudydirs{sim},'solve'),'dir')
      delete(fullfile(mystudydirs{sim},'solve','core*'));
      delete(fullfile(mystudydirs{sim},'solve','params.mat'));
      delete(fullfile(mystudydirs{sim},'solve','solve_ode*'));
      %delete(fullfile(mystudydirs{sim},'solve','*'))
      %rmdir(fullfile(mystudydirs{sim},'solve'));
      %rmdir(fullfile(mystudydirs{sim}));
    end
  end
  
  % Try to force removal of parfor folders (this might produce error if
  % lock files are present
  try rmdir(fullfile(options.study_dir,'output_parfor*'),'s');
  catch
      %warning('If this is reached, it is likely because there are locked files in the folder, which cannot be deleted. The error is /.nfs00000000a617655e0000c15a: Device or resource busy. Fix for this is not obvious');
  end
  
  % Delete any core files in parent directory
  delete(fullfile(options.study_dir,'core*'));
  
  % Verify all core files are deleted
  [~,result] = system('find * -name "core*"','-echo');
  if ~isempty(result); fprintf(strcat(result,'\n')); warning('Core files found. Consider deleting to free up space'); end
  
  % close pool
%   delete(gcp)
% todo: sort data sets by things varied in modifications_set
% todo: Figure out how to delete locked .nfs files
  return
end

% 1.4 prepare study_dir and studyinfo if saving data
if isempty(options.studyinfo)
  [studyinfo,options]=SetupStudy(model,'modifications_set',modifications_set,'simulator_options',options,'process_id',options.sim_id);
else
  studyinfo=options.studyinfo;
end

% put solver blocks in try statement to catch and handle errors/cleanup
cwd=pwd; % record current directory
try
  
  % functions:          outputs:
  % WriteDynaSimSolver    m-file for DynaSim solver
  % WriteMatlabSolver   m-file for Matlab solver (including @odefun)
  % PrepareMEX          mex-file for m-file

  % 1.5 loop over simulations, possibly varying things
  base_model=model;
  data_index=0;
  for sim=1:length(modifications_set)
    if ~strcmp(pwd,cwd) % move back to original directory before potentially regenerating to make sure the model files used are the same
      cd(cwd);
    end
    % get index for this simulation
    if ~isempty(options.sim_id)
      sim_ind=find([studyinfo.simulations.sim_id]==options.sim_id);
      sim_id=options.sim_id;
    else
      sim_ind=sim;
      sim_id=sim;
    end
    if options.save_data_flag
      % check if output data already exists. load if so and skip simulation
      data_file=studyinfo.simulations(sim_ind).data_file;
      if exist(data_file,'file') && options.overwrite_flag==0
        if 1%options.verbose_flag
          % note: this is important, should always display
          fprintf('loading data from %s\n',data_file);
        end
        tmpdata=ImportData(data_file,'process_id',sim_id);
        update_data; % concatenate data structures across simulations
        continue; % skip to next simulation
      end
    end
    % apply modifications for this point in search space
    if ~isempty(modifications_set{sim})
      model=ApplyModifications(base_model,modifications_set{sim});
    end
    % update studyinfo
    if options.save_data_flag
      %studyinfo=UpdateStudy(studyinfo.study_dir,'process_id',sim_id,'status','started','model',model,'simulator_options',options,'verbose_flag',options.verbose_flag);
    end
    % execute experiment
    if isa(options.experiment,'function_handle')
      % EXPERIMENT (wrapping around a set of simulations)
      if options.cluster_flag && options.compile_flag
        warning('compiled solver is not available for experiments on the cluster. Simulation will be run in Matlab.');
      end
      % from varargin...
      % remove 'experiment', 'modifications', 'vary', 'cluster_flag' to avoid undesired recursive action in experiment function
      % remove 'save_data_flag' to prevent individual simulations from being saved during experiment
      keyvals=RemoveKeyval(varargin,{'experiment','cluster_flag','vary','modifications','save_data_flag'});
      if ~isempty(options.experiment_options)
        % user-supplied experiment options override any found in SimulateModel options
        keyvals=RemoveKeyval(keyvals,options.experiment_options(1:2:end));
        keyvals=cat(2,keyvals,options.experiment_options);
      end
      tmpdata=feval(options.experiment,model,keyvals{:});
    else
      % NOT AN EXPERIMENT (single simulation)
      %% 2.0 prepare solver function (solve_ode.m/mex)
      % - matlab solver: create @odefun with vectorized state variables
      % - DynaSim solver: write solve_ode.m and params.mat  (based on dnsimulator())    
      % check if model solver needs to be created 
      % (i.e., if is first simulation or a search space varying mechanism list)
      if sim==1 || (~isempty(modifications_set{1}) && any(cellfun(@(x)strcmp(x{2},'mechanism_list'),modifications_set)))
        % prepare file that solves the model system
        if isempty(options.solve_file) || (~exist(options.solve_file,'file') && ~exist([options.solve_file '.mexa64'],'file') &&  ~exist([options.solve_file '.mexa32'],'file') && ~exist([options.solve_file '.mexmaci64'],'file'))
          options.solve_file=GetSolveFile(model,studyinfo,options); % store name of solver file in options struct
        end
        % todo: consider providing better support for studies that produce different m-files per sim (e.g., varying mechanism_list)
        if options.verbose_flag
          fprintf('SIMULATING MODEL:\n');
          fprintf('solving system using %s\n',options.solve_file);
        end
      else
        % use previous solve_file
      end    
      [fpath,fname,fext]=fileparts(options.solve_file);

      %% 3.0 integrate model with solver of choice and prepare output data
      % - matlab solver: solve @odefun with feval and solver_options
      % - DynaSim solver: run solve_ode.m or create/run MEX
      % move to directory with solver file
      if options.verbose_flag
        fprintf('changing directory to %s\n',fpath);
      end
      cd(fpath);
      % save parameters there
      warning('off','catstruct:DuplicatesFound');
      p=catstruct(CheckSolverOptions(options),model.parameters);
      param_file=fullfile(fpath,'params.mat');
      if options.verbose_flag
        fprintf('saving model parameters: %s\n',param_file);
      end
      %pause(.01);
      % solve system
      if options.disk_flag  % ### data stored on disk during simulation ###
        sim_start_time=tic;
        save(param_file,'p'); % save params immediately before solving
        csv_data_file=feval(fname);  % returns name of file storing the simulated data
        duration=toc(sim_start_time);
        if nargout>0 || options.save_data_flag
          tmpdata=ImportData(csv_data_file,'process_id',sim_id); % eg, data.csv
        end
      else                  % ### data stored in memory during simulation ###
        % create list of output variables to capture
        output_variables=cat(2,'time',model.state_variables);
        if ~isempty(model.monitors)
          output_variables=cat(2,output_variables,fieldnames(model.monitors)');
        end
        if ~isempty(model.fixed_variables)
          fields=fieldnames(model.fixed_variables)';
          output_variables=cat(2,output_variables,fields);
          num_fixed_variables=length(fields);
        else
          num_fixed_variables=0;
        end
        % run simulation
        if options.verbose_flag
          fprintf('Running simulation %g/%g (solver=''%s'', dt=%g, tspan=[%g %g]) ...\n',sim,length(modifications_set),options.solver,options.dt,options.tspan);
        end
        sim_start_time=tic;
        outputs=cell(1,length(output_variables)); % preallocate for PCT compatibility
        save(param_file,'p'); % save params immediately before solving
        [outputs{1:length(output_variables)}]=feval(fname);
        duration=toc(sim_start_time);
        % prepare DynaSim data structure
        % organize simulated data in data structure (move time to last)
        tmpdata.labels=output_variables([2:length(output_variables)-num_fixed_variables 1]);
        for i=1:length(output_variables)
          if ~isempty(model.fixed_variables) && isfield(model.fixed_variables,output_variables{i})
            % store fixed variables in model substructure
            model.fixed_variables.(output_variables{i})=outputs{i};
          else
            % store state variables and monitors as data fields
            tmpdata.(output_variables{i})=outputs{i};
          end
          outputs{i}=[]; % clear assigned outputs from memory
        end
      end
      if options.verbose_flag
        fprintf('Elapsed time: %g seconds.\n',duration); 
      end
      % add metadata to tmpdata
      tmpdata.simulator_options=options; % store simulator controls
      if options.store_model_flag==1  % optionally store the simulated model
        tmpdata.model=model;
      end
    end
    prepare_varied_metadata;
    % save single data set and update studyinfo
    if options.save_data_flag
      ExportData(tmpdata,'filename',data_file,'format','mat','verbose_flag',options.verbose_flag);
      %studyinfo=UpdateStudy(studyinfo.study_dir,'process_id',sim_id,'status','finished','duration',duration,'solve_file',options.solve_file,'email',options.email,'verbose_flag',options.verbose_flag,'model',model,'simulator_options',options);
    end
    % do post-simulation analysis and plotting
    if ~isempty(options.analysis_functions) || ~isempty(options.plot_functions)
      if options.save_data_flag || options.save_results_flag
        % do analysis and plotting while saving results
        siminfo=studyinfo.simulations(sim_ind);
        for f=1:length(siminfo.result_functions)
          result=AnalyzeData(tmpdata,siminfo.result_functions{f},'result_file',siminfo.result_files{f},'save_data_flag',1,'save_results_flag',1,siminfo.result_options{f}{:});
          % since the plots are saved, close all generated figures
          if all(ishandle(result))
            close(result);
          end
        end
      else
        % do analysis and plotting without saving results
        if ~isempty(options.analysis_functions)
          for f=1:length(options.analysis_functions)
            tmpdata=AnalyzeData(tmpdata,options.analysis_functions{f},'result_file',[],'save_data_flag',0,'save_results_flag',options.save_results_flag,options.analysis_options{f}{:});
          end
        end
        if ~isempty(options.plot_functions)
          for f=1:length(options.plot_functions)
            AnalyzeData(tmpdata,options.plot_functions{f},'result_file',[],'save_data_flag',0,'save_results_flag',options.save_results_flag,options.plot_options{f}{:});
          end
        end
      end
    end
    if nargout>0
      update_data; % concatenate data structures across simulations
    end
  end % end loop over sims
  cleanup('success');
catch err % error handling
  if options.compile_flag && ~isempty(options.solve_file)
    if options.verbose_flag
      fprintf('removing failed compiled solve file: %s\n',options.solve_file);
    end
    delete([options.solve_file '*']);
  end
  DisplayError(err);
  %keyboard
  % update studyinfo
  if options.save_data_flag
    studyinfo=UpdateStudy(studyinfo.study_dir,'process_id',sim_id,'status','failed','verbose_flag',options.verbose_flag);
    data=studyinfo;
  end
  cleanup('error');
  return  
end

% ---------------------------------------------
% todo:
% - create function that constructs @odefun
% - add support for built-in matlab solvers
% - create helper function that handles log files (creation, standardized format,...)
% ---------------------------------------------

% -------------------------
% NESTED FUNCTIONS
% -------------------------
  function update_data
    % store tmpdata
    if sim==1
      % replicate first data set as preallocation for all
      data=repmat(tmpdata,[1 length(modifications_set)]);
      data_index=length(tmpdata);
    else
      inds=data_index+(1:length(tmpdata)); % support multiple data sets returned by experiments
      data(inds)=tmpdata;
      data_index=inds(end);
    end
  end
  function prepare_varied_metadata
    % add things varied to tmpdata
    mods={};
    if ~isempty(options.modifications)
      mods=cat(1,mods,expand_modifications(options.modifications));
    end
    if ~isempty(modifications_set{sim})
      tmp_mods=expand_modifications(modifications_set{sim});
      mods=cat(1,mods,tmp_mods);
    end
    if isa(options.experiment,'function_handle')
      for j=1:length(tmpdata)
        tmpdata(j).simulator_options.modifications=mods;
      end
    end
    if ~isempty(mods)
      if isfield(tmpdata,'varied')
        varied=tmpdata(1).varied;
      else
        varied={};
      end
      for ii=1:size(mods,1)
        % prepare valid field name for thing varied:
        fld=[mods{ii,1} '_' mods{ii,2}];
        % convert arrows and periods to underscores
        fld=regexprep(fld,'(->)|(<-)|(-)|(\.)','_');
        % remove brackets and parentheses
        fld=regexprep(fld,'[\[\]\(\)\{\}]','');
        for j=1:length(tmpdata)
          tmpdata(j).(fld)=mods{ii,3};
        end
        if ~ismember(fld,varied)
          varied{end+1}=fld;
        end
      end
      for j=1:length(tmpdata)
        tmpdata(j).varied=varied;
      end
    end    
    % convert tmpdata to single precision
    if strcmp(options.precision,'single')
      for j=1:length(tmpdata)
        for k=1:length(tmpdata(j).labels)
          fld=tmpdata(j).labels{k};
          tmpdata(j).(fld)=single(tmpdata(j).(fld));
        end
      end    
    end
  end
    
  function cleanup(status)
      % remove temporary files and optionally store info for debugging
      % ...
      % return to original directory
      if options.verbose_flag
        fprintf('changing directory to %s\n',cwd);
      end
      cd(cwd);
    switch status
      case 'success'
        % ...
      case 'error'
        % ... error logs
    end
  end

  function all_ICs=ProcessNumericICs
    % first, figure out how many IC values we need (i.e., how many state
    % variables we need across all cells).    
    var_names=model.state_variables;
    [nvals_per_var,monitor_counts]=GetOutputCounts(model);
    num_state_variables=sum(nvals_per_var);
    % check that the correct number of IC values was provided
    if length(options.ic)~=num_state_variables
      error('incorrect number of initial conditions. %g values are needed for %g state variables across %g cells',num_state_variables,length(model.state_variables),sum(pop_sizes));
    end
    % organize user-supplied ICs into array for each state variable (assume
    cnt=0; all_ICs=[];
    for i=1:length(var_names)
      ICs=options.ic(cnt+(1:nvals_per_var(i)));
      % store ICs as string for writing solve_ode and consistent evaluation
      all_ICs.(var_names{i})=sprintf('[%s]',num2str(ICs));
      cnt=cnt+nvals_per_var(i);
    end    
  end
    
end

function modifications=expand_modifications(mods)
  % purpose: expand simultaneous modifications into larger list
  modifications={};
  for i=1:size(mods,1)
    % get object list without grouping symbols: ()[]{}
    objects=regexp(mods{i,1},'[^\(\)\[\]\{\},]+','match');
    variables=regexp(mods{i,2},'[^\(\)\[\]\{\},]+','match');
    for j=1:length(objects)
      for k=1:length(variables)
        modifications(end+1,1:3)={objects{j},variables{k},mods{i,3}};
      end
    end
  end
end

function [model,options]=extract_vary_statement(model,options)
  % purpose: extract vary statement, remove from model, and set options.vary
  if ischar(model) && any(regexp(model,';\s*vary\(.*\)','once'))
    % extract vary statement
    str=regexp(model,';\s*(vary\(.*\);?)','tokens','once');
    % remove from model
    model=strrep(model,str{1},'');
    % set options
    var=regexp(str{1},'\((.*)=','tokens','once'); % variable
    val=regexp(str{1},'=(.*)\)','tokens','once'); % values
    options.vary={'pop1',var{1},eval(val{1})};
  end
end

function options = backward_compatibility(options)
% option_names: (old_name, new_name; ...}
option_names = {...
  'override','modifications';
  'timelimits','tspan';
  'IC','ic';
  'verbose','verbose_flag';
  'SOLVER','solver';
  'nofunctions','reduce_function_calls_flag';
  'dsfact','downsample_factor';
  'memlimit','memory_limit';
  };
if any(ismember(option_names(:,1),options(1:2:end)))
  for i=1:size(option_names,1)
    % check if any options have this old name
    if ismember(option_names{i,1},options(1:2:end))
      ind=find(ismember(options(1:2:end),option_names{i,1}));
      % replace old option name by new option name
      options{2*ind-1}=option_names{i,2};
    end
  end
end
end

