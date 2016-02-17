function [outfile,options]=WriteDynaSimSolver(model,varargin)
%% solver_file=WriteDynaSimSolver(model,varargin)
% Purpose: write m-file that numerically inteegrates the model
% Inputs:
%   model: DynaSim model structure (see GenerateModel)
%   options:
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
%     'disk_flag'     : whether to write to disk during simulation instead of storing in memory {0 or 1} (default: 0)
% 
% Outputs:
%   solver_file (e.g., solve_ode.m): function that numerically integrates the model
% 
% Example 1: test solver production, display function in standard output 
% model=GenerateModel; % test model
% without writing anything to disk:
% WriteDynaSimSolver(model,'fileID',1,'save_parameters_flag',0,'solver','rk4');
% WriteDynaSimSolver(model,'fileID',1,'save_parameters_flag',1,'solver','rk4');
% model=PropagateFunctions(model);
% WriteDynaSimSolver(model,'fileID',1,'save_parameters_flag',0,'solver','rk4');
% WriteDynaSimSolver(model,'fileID',1,'save_parameters_flag',1,'solver','rk4');
% 
% Example 2: real-time downsampling
% WriteDynaSimSolver(model,'downsample_factor',10,'fileID',1,'solver','rk4');
% 
% Example 3: real-time writing data to disk (reduce memory demand)
% WriteDynaSimSolver(model,'disk_flag',1,'fileID',1,'solver','rk4');
% 
% Example 4: maximally conserve memory: downsample and write to disk
% WriteDynaSimSolver(model,'disk_flag',1,'downsample_factor',10,'fileID',1,'solver','rk4');
% 
% Examples:
% WriteDynaSimSolver(model,'solver','euler');
% WriteDynaSimSolver(model,'solver','rk2');
% WriteDynaSimSolver(model,'solver','rk4');
% 
% model=GenerateModel; % test model
% tic; WriteDynaSimSolver(model,'save_parameters_flag',0,'reduce_function_calls_flag',0,'solver','rk4','filename','solve_ode.m'); v=solve_ode; plot(v); toc
% tic; WriteDynaSimSolver(model,'save_parameters_flag',1,'reduce_function_calls_flag',0,'solver','rk4','filename','solve_ode.m'); v=solve_ode; plot(v); toc
% tic; WriteDynaSimSolver(model,'save_parameters_flag',0,'reduce_function_calls_flag',1,'solver','rk4','filename','solve_ode.m'); v=solve_ode; plot(v); toc
% tic; WriteDynaSimSolver(model,'save_parameters_flag',1,'reduce_function_calls_flag',1,'solver','rk4','filename','solve_ode.m'); v=solve_ode; plot(v); toc
% tic; WriteDynaSimSolver(model,'save_parameters_flag',1,'reduce_function_calls_flag',1,'solver','rk2','filename','solve_ode.m'); v=solve_ode; plot(v); toc
% tic; WriteDynaSimSolver(model,'save_parameters_flag',1,'reduce_function_calls_flag',1,'solver','rk1','filename','solve_ode.m'); v=solve_ode; plot(v); toc
% tic; WriteDynaSimSolver(model,'save_parameters_flag',1,'reduce_function_calls_flag',1,'solver','rk1','filename','solve_ode.m','downsample_factor',10); v=solve_ode; plot(v); toc
% tic; WriteDynaSimSolver(model,'save_parameters_flag',1,'reduce_function_calls_flag',1,'solver','rk2','filename','solve_ode.m','dt',.001,'downsample_factor',10); v=solve_ode; plot(v); toc
% tic; WriteDynaSimSolver(model,'save_parameters_flag',1,'reduce_function_calls_flag',1,'solver','rk1','filename','solve_ode.m','dt',.001,'downsample_factor',10); v=solve_ode; plot(v); toc
% tic; WriteDynaSimSolver(model,'save_parameters_flag',1,'reduce_function_calls_flag',1,'solver','rk1','filename','solve_ode.m','dt',.005,'tspan',[0 200],'downsample_factor',10); v=solve_ode; plot(v); toc
% 
% See also: GetSolveFile, SimulateModel, WriteMatlabSolver

% dependencies: CheckOptions, CheckModel

% Check inputs
options=CheckOptions(varargin,{...
  'tspan',[0 100],[],...          % [beg,end] (units must be consistent with dt and equations)
  'downsample_factor',1,[],...    % downsampling applied during simulation (only every downsample_factor-time point is stored in memory or written to disk)
  'random_seed','shuffle',[],...        % seed for random number generator (usage: rng(random_seed))
  'solver','rk4',{'euler','rk1','rk2','rk4','modified_euler','rungekutta','rk'},... % DynaSim and built-in Matlab solvers
  'disk_flag',0,{0,1},...            % whether to write to disk during simulation instead of storing in memory
  'dt',.01,[],...                 % time step used for fixed step DynaSim solvers
  'reduce_function_calls_flag',1,{0,1},...   % whether to eliminate internal (anonymous) function calls
  'save_parameters_flag',1,{0,1},...
  'filename',[],[],...         % name of solver file that integrates model
  'data_file','data.csv',[],... % name of data file if disk_flag=1
  'fileID',1,[],...
  'compile_flag',exist('codegen')==6,{0,1},... % whether to prepare script for being compiled using coder instead of interpreting Matlab
  'verbose_flag',1,{0,1},...
  },false);
model=CheckModel(model); 
separator=','; % ',', '\\t'

%% 1.0 prepare model info
parameter_prefix='p.';%'pset.p.';
state_variables=model.state_variables;

% 1.1 eliminate internal (anonymous) function calls from model equations
if options.reduce_function_calls_flag==1
  model=PropagateFunctions(model);
end

% 1.1 prepare parameters
if options.save_parameters_flag
  % add parameter struct prefix to parameters in model equations
  model=PropagateParameters(model,'action','prepend','prefix',parameter_prefix);
  % set and capture numeric seed value
  if options.compile_flag==1
    % todo: make seed string (eg, 'shuffle') from param struct work with coder (options.compile_flag=1)
    % (currently raises error: "String input must be constant")
    % workaround: (shuffle here and get numeric seed for MEX-compatible params.mat)
    rng(options.random_seed);
    options.random_seed=getfield(rng,'Seed');  % <-- current active seed
  end
  % set parameter file name (save with m-file)
  [fpath,fname,fext]=fileparts(options.filename);
  param_file_name = fullfile(fpath,'params.mat');
  % save parameters to disk
  warning('off','catstruct:DuplicatesFound');
  p=catstruct(CheckSolverOptions(options),model.parameters);
  if options.verbose_flag
    fprintf('saving params.mat\n');
  end
  save(param_file_name,'p');
else
  % insert parameter values into model expressions
  model=PropagateParameters(model,'action','substitute');
end

% 1.2 prepare list of outputs (state variables and monitors)
tmp=cellfun(@(x)[x ','],model.state_variables,'uni',0);
tmp=[tmp{:}];
output_string=tmp(1:end-1);
if ~isempty(model.monitors)
  tmp=cellfun(@(x)[x ','],fieldnames(model.monitors),'uni',0);
  tmp=[tmp{:}];
  output_string=[output_string ',' tmp(1:end-1)];
end
if ~isempty(model.fixed_variables)
  tmp=cellfun(@(x)[x ','],fieldnames(model.fixed_variables),'uni',0);
  tmp=[tmp{:}];
  output_string=[output_string ',' tmp(1:end-1)];
end  
output_string=['[T,' output_string ']']; % state vars, monitors, time vector

%% 2.0 write m-file solver
% 2.1 create mfile
if ~isempty(options.filename)
  if options.verbose_flag
    fprintf('creating solver file: %s\n',options.filename);
  end
  fid=fopen(options.filename,'wt');
else
  fid=options.fileID;
end
outfile=fopen(fid);

if options.disk_flag==1
  fprintf(fid,'function data_file=solve_ode\n');
  % create output data file
  fprintf(fid,'%% ------------------------------------------------------------\n');
  fprintf(fid,'%% Open output data file:\n');
  fprintf(fid,'data_file=''%s'';\n',options.data_file);
  %fprintf(fid,'fileID=fopen(data_file,''wt'');\n'); % <-- 'wt' does not
  % compile in linux. 't' may be necessary on PC. may need to look into this
  fprintf(fid,'fileID=fopen(data_file,''w'');\n');
  % write headers
  [state_var_counts,monitor_counts]=GetOutputCounts(model);
  fprintf(fid,'fprintf(fileID,''time%s'');\n',separator);
  if ~isempty(model.state_variables)
    for i=1:length(model.state_variables)
      fprintf(fid,'for i=1:%g, fprintf(fileID,''%s%s''); end\n',state_var_counts(i),model.state_variables{i},separator);    
    end
  end
  if ~isempty(model.monitors)
    monitor_names=fieldnames(model.monitors);
    for i=1:length(monitor_names)
      fprintf(fid,'for i=1:%g, fprintf(fileID,''%s%s''); end\n',monitor_counts(i),monitor_names{i},separator);    
    end
  end
  fprintf(fid,'fprintf(fileID,''\\n'');\n');
else
  fprintf(fid,'function %s=solve_ode\n',output_string);
end

% 2.3 load parameters
if options.save_parameters_flag
  fprintf(fid,'%% ------------------------------------------------------------\n');
  fprintf(fid,'%% Parameters:\n');
  fprintf(fid,'%% ------------------------------------------------------------\n');
  fprintf(fid,'p=load(''params.mat'',''p''); p=p.p;\n');
end

% write tspan, dt, and downsample_factor
if options.save_parameters_flag
  fprintf(fid,'downsample_factor=%sdownsample_factor;\n',parameter_prefix);
  fprintf(fid,'dt=%sdt;\n',parameter_prefix);
  fprintf(fid,'T=(%stspan(1):dt:%stspan(2))'';\n',parameter_prefix,parameter_prefix);
else
  fprintf(fid,'tspan=[%g %g];\ndt=%g;\ndownsample_factor=%g;\n',options.tspan,options.dt,options.downsample_factor);
  fprintf(fid,'T=(tspan(1):dt:tspan(2))'';\n');
end
% write calculation of time vector and derived parameters: ntime, nsamp
fprintf(fid,'ntime=length(T);\nnsamp=length(1:downsample_factor:ntime);\n');

% 2.4 evaluate fixed variables
if ~isempty(model.fixed_variables)
  fprintf(fid,'%% ------------------------------------------------------------\n');
  fprintf(fid,'%% Fixed variables:\n');
  fprintf(fid,'%% ------------------------------------------------------------\n');
  names=fieldnames(model.fixed_variables);
  expressions=struct2cell(model.fixed_variables);
  for i=1:length(names)
    fprintf(fid,'%s = %s;\n',names{i},expressions{i});
  end
end

% 2.5 evaluate function handles
if ~isempty(model.functions) && options.reduce_function_calls_flag==0
  fprintf(fid,'%% ------------------------------------------------------------\n');
  fprintf(fid,'%% Functions:\n');
  fprintf(fid,'%% ------------------------------------------------------------\n');
  names=fieldnames(model.functions);
  expressions=struct2cell(model.functions);
  for i=1:length(names)
    fprintf(fid,'%s = %s;\n',names{i},expressions{i});
  end
end

% 2.6 prepare storage
fprintf(fid,'%% ------------------------------------------------------------\n');
fprintf(fid,'%% Initial conditions:\n');
fprintf(fid,'%% ------------------------------------------------------------\n');
% 2.2 set random seed
fprintf(fid,'%% seed the random number generator\n');
if options.save_parameters_flag
  fprintf(fid,'rng(%srandom_seed);\n',parameter_prefix);
else  
  if ischar(options.random_seed)
    fprintf(fid,'rng(''%s'');\n',options.random_seed);
  elseif isnumeric(options.random_seed)
    fprintf(fid,'rng(%g);\n',options.random_seed);
  end
end
% initialize time
fprintf(fid,'t=0; k=1;\n');

% todo: get coder varsize working with new format:

% prepare for compilation
% if options.compile_flag==1
%   for i = 1:length(state_variables)
%     %fprintf(fid,'coder.varsize(''%s'',[1e4,1],[true,false]);\n',state_variables{i}); % population size up to 1e4 for variable initial conditions
%     fprintf(fid,'coder.varsize(''%s'',[1e8,1e4],[true,true]);\n',state_variables{i}); % population size up to 1e4 for variable initial conditions
%   end  
% end

% preallocate and initialize matrices to store results
if options.disk_flag==1
  % add time zero to output data file
  fprintf(fid,'fprintf(fileID,''0%s'');\n',separator);
end
fprintf(fid,'%% STATE_VARIABLES:\n');
% STATE_VARIABLES
nvals_per_var=zeros(1,length(state_variables));
IC_expressions=struct2cell(model.ICs);
for i=1:length(state_variables)
  % initialize var_last
  if options.downsample_factor>1 || options.disk_flag==1
    % set var_last=IC;
    fprintf(fid,'%s_last = %s;\n',state_variables{i},IC_expressions{i});
  end
  if options.disk_flag==1
    % print var_last
    var_last=sprintf('%s_last',state_variables{i});
    fprintf(fid,'for i=1:numel(%s), fprintf(fileID,''%%g%s'',%s(i)); end\n',var_last,separator,var_last);
  else
    % preallocate state variables
    parts=regexp(state_variables{i},'_','split');
    if numel(parts)==4 % has connection mechanism namespace: target_source_mechanism
      % state variables defined in connection mechanisms are assumed to
      % have dimensionality of the source population
      part=parts{2};
    else % has intrinsic mechanism or population namespace: target_mechanism
      % state variables defined in intrinsic mechanisms or population
      % equations have dimensionality of the target population
      part=parts{1};
    end
    if options.save_parameters_flag
      fprintf(fid,'%s = zeros(nsamp,%s%s_Npop);\n',state_variables{i},parameter_prefix,part);
    else
      fprintf(fid,'%s = zeros(nsamp,%g);\n',state_variables{i},model.parameters.([part '_Npop']));
    end
    nvals_per_var(i)=model.parameters.([part '_Npop']);
    % initialize state variables
    if options.downsample_factor==1
      % set var(1,:)=IC;
      fprintf(fid,'  %s(1,:) = %s;\n',state_variables{i},IC_expressions{i});      
    else
      % set var(1,:)=var_last;
      fprintf(fid,'  %s(1,:) = %s_last;\n',state_variables{i},state_variables{i});      
    end
  end
end
if ~isempty(model.monitors)
  fprintf(fid,'%% MONITORS:\n');
  % MONITORS    
  monitor_names=fieldnames(model.monitors);
  monitor_expression=struct2cell(model.monitors);
  for i=1:length(monitor_names)
    if ~isempty(regexp(monitor_names{i},'_spikes$','once'))
    % set expression if monitoring spikes
      if options.disk_flag==1
        error('spike monitoring is not supported for writing data to disk at this time.');
        % todo: add support for spike monitoring with data written to
        % disk. approach: load data at end of simulation and do post-hoc
        % spike finding. if saving *.mat, use matfile() to load only the
        % variables in which to search. add spikes to output data file.
      end
      if isempty(monitor_expression{i})
        spike_threshold=0;
      elseif isempty(regexp(monitor_expression{i},'[^\d]','once'))
        spike_threshold=str2num(monitor_expression{i});
        monitor_expression{i}=[];
      else
        spike_threshold=monitor_expression{i};
        monitor_expression{i}=[];
      end
      % approach: add conditional check for upward threshold crossing
      var_spikes=regexp(monitor_names{i},'(.*)_spikes$','tokens','once');
      var_spikes=var_spikes{1}; % variable to monitor
      if ismember(var_spikes,model.state_variables)
        model.conditionals(end+1).namespace='spike_monitor';
        if isnumeric(spike_threshold)
          model.conditionals(end).condition=...
            sprintf('%s(n,:)>=%g&%s(n-1,:)<%g',var_spikes,spike_threshold,var_spikes,spike_threshold);
        else
          model.conditionals(end).condition=...
            sprintf('%s(n,:)>=%s&%s(n-1,:)<%s',var_spikes,spike_threshold,var_spikes,spike_threshold);
        end
        model.conditionals(end).action=sprintf('%s(n,conditional_test)=1',monitor_names{i});
        model.conditionals(end).else=[];
        % move spike monitor to first position (ie.., to evaluate before other conditionals)
        model.conditionals=model.conditionals([length(model.conditionals) 1:length(model.conditionals)-1]);
        % remove from monitor list
        model.monitors=rmfield(model.monitors,monitor_names{i});
      end      
    elseif isempty(monitor_expression{i}) && isfield(model.functions,monitor_names{i})
    % set expression if monitoring function referenced by name
      tmp=regexp(model.functions.(monitor_names{i}),'@\([a-zA-Z][\w,]*\)\s*(.*)','tokens','once');
      monitor_expression{i}=tmp{1};
      model.monitors.(monitor_names{i})=tmp{1};
    end
    % initialize mon_last
    if options.downsample_factor>1 || options.disk_flag==1
      % set mon_last=f(IC);
      tmp=cell2struct({monitor_expression{i}},{monitor_names{i}},1);
      print_monitor_update(fid,tmp,'_last',state_variables,'_last');
    end
    if options.disk_flag==1
      % print mon_last
      mon_last=sprintf('%s_last',monitor_names{i});
      fprintf(fid,'for i=1:numel(%s), fprintf(fileID,''%%g%s'',%s(i)); end\n',mon_last,separator,mon_last);
    else
      % preallocate monitors
      parts=regexp(monitor_names{i},'_','split');
      if options.save_parameters_flag
        fprintf(fid,'%s = zeros(nsamp,%s%s_Npop);\n',monitor_names{i},parameter_prefix,parts{1});
      else
        fprintf(fid,'%s = zeros(nsamp,%g);\n',monitor_names{i},model.parameters.([parts{1} '_Npop']));
      end    
      if isempty(monitor_expression{i})
        continue;
      end
      % initialize monitors
      if options.downsample_factor==1
        % set mon(1,:)=f(IC);
        tmp=cell2struct({monitor_expression{i}},{monitor_names{i}},1);
        print_monitor_update(fid,tmp,'(1,:)',state_variables,'(1,:)');
      else
        % set mon(1,:)=mon_last;
        tmp=cell2struct({monitor_expression{i}},{monitor_names{i}},1);
        print_monitor_update(fid,tmp,'(1,:)',state_variables,'_last');
      end
    end
  end
end
if options.disk_flag==1
  % go to new line for next time point
  fprintf(fid,'fprintf(fileID,''\\n'');\n');
end

% determine how to index each state variable based on how often state
% variables are stored, whether they are written to disk or stored in
% memory, and whether the state variable matrix is one- or two-dimensional
index_lasts=cell(1,length(state_variables));
index_nexts=cell(1,length(state_variables));
index_temps=repmat({'_last'},[1 length(state_variables)]);
for i=1:length(state_variables)
  if options.downsample_factor==1 && options.disk_flag==0
    % store state directly into state variables on each integration step
    if nvals_per_var(i)>1 % use full 2D matrix indexing
      index_lasts{i}='(n-1,:)';
      index_nexts{i}='(n,:)';
    else % use more concise 1D indexing because it is much faster for some Matlab-specific reason...
      index_lasts{i}='(n-1)';
      index_nexts{i}='(n)';
    end  
  elseif options.downsample_factor>1 && options.disk_flag==0
    % store state in var_last then update state variables on each downsample_factor integration step
    index_lasts{i}='_last';
    if nvals_per_var(i)>1
      index_nexts{i}='(n,:)';
    else
      index_nexts{i}='(n)';
    end
  elseif options.disk_flag==1
    % always store state in var_last and write on each downsample_factor integration step
      index_lasts{i}='_last';
      index_nexts{i}='_last';  
  end
end

% add index to state variables in ODEs
odes = struct2cell(model.ODEs);
for i=1:length(odes)
  for j=1:length(state_variables)
    odes{i}=dynasim_strrep(odes{i},state_variables{j},[state_variables{j} index_lasts{j}]);
  end
end

% write code to do numerical integration
fprintf(fid,'%% ###########################################################\n');
fprintf(fid,'%% Numerical integration:\n');
fprintf(fid,'%% ###########################################################\n');
fprintf(fid,'n=2;\n');
fprintf(fid,'for k=2:ntime\n'); % time index
fprintf(fid,'  t=T(k-1);\n');
if options.downsample_factor==1 && options.disk_flag==0 % store every time point, do not use var_last
  % update_vars;      % var(:,k-1)->var(k,:) or var(k-1)->var(k)
  update_vars(index_nexts);
  % conditionals;     % var(k,:)->var(k,:) or var(k)->var(k)
  print_conditional_update(fid,model.conditionals,index_nexts,state_variables)
  % update_monitors;  % mon(:,k-1)->mon(k,:) or mon(k-1)->mon(k)
  print_monitor_update(fid,model.monitors,index_nexts,state_variables);
  fprintf(fid,'  n=n+1;\n');
else % store every downsample_factor time point in memory or on disk
  % update_vars;      % var_last->var_last
  update_vars(index_temps);
  % conditionals;     % var_last->var_last
  print_conditional_update(fid,model.conditionals,index_temps,state_variables);
  % check if it is time to store data
  fprintf(fid,'if mod(k,downsample_factor)==0 %% store this time point\n');
  if options.disk_flag==1 % write to disk
    fprintf(fid,'  %% ------------------------------------------------------------\n');
    fprintf(fid,'  %% Write state variables and monitors to disk:\n');
    fprintf(fid,'  %% ------------------------------------------------------------\n');
    % print current time point to data file
    fprintf(fid,'  fprintf(fileID,''%%g%s'',T(k));\n',separator);  
    % print state variables
    for i=1:length(state_variables)
      var_last=sprintf('%s_last',state_variables{i});
      fprintf(fid,'  for i=1:numel(%s), fprintf(fileID,''%%g%s'',%s(i)); end\n',var_last,separator,var_last);
    end
    % print monitors
    if ~isempty(model.monitors)
      for i=1:length(monitor_names)
          mon_last=sprintf('%s_last',monitor_names{i});
          fprintf(fid,'  for i=1:numel(%s), fprintf(fileID,''%%g%s'',%s(i)); end\n',mon_last,separator,mon_last);
      end
    end
    fprintf(fid,'  fprintf(fileID,''\\n'');\n');
  else            % store in memory
    % update_vars;    % var_last -> var(n,:) or var(n)
    print_var_update_last(fid,index_nexts,state_variables)
    % update_monitors;% f(var_last) -> mon(n,:) or mon(n)
    print_monitor_update(fid,model.monitors,index_nexts,state_variables);
  end
  fprintf(fid,'  n=n+1;\n');
  fprintf(fid,'end\n');
end
fprintf(fid,'end\n');
fprintf(fid,'T=T(1:downsample_factor:ntime);\n');

% cleanup
if options.disk_flag==1
  fprintf(fid,'%% ------------------------------------------------------------\n');
  fprintf(fid,'%% Close output data file:\n');
  fprintf(fid,'fclose(fileID);\n');
  fprintf(fid,'%% ------------------------------------------------------------\n');
end
if ~strcmp(outfile,'"stdout"')
  fclose(fid);
  % wait for file before continuing to simulation
  while ~exist(outfile,'file')
    pause(.01);
  end
end

  % ########################################
  % NESTED FUNCTIONS
  % ########################################
  function update_vars(index_nexts_)
    switch options.solver
      case {'euler','rk1'}
        print_k(fid,odes,'_k1',state_variables);                              % write k1 using model.ODEs
        % write update for state variables
        print_var_update(fid,index_nexts_,index_lasts,...
          'dt*%s_k1',state_variables);
      case {'rk2','modified_euler'}
        print_k(fid,odes,'_k1',state_variables);                              % write k1 using model.ODEs
        odes_k2=update_odes(odes,'_k1','.5*dt',state_variables,index_lasts);   % F(*,yn+dt*k1/2)
        fprintf(fid,'  t=t+.5*dt;\n');                                          % update t before writing k2
        print_k(fid,odes_k2,'_k2',state_variables);                           % write k2 using odes_k2
        % write update for state variables
        print_var_update(fid,index_nexts_,index_lasts,...
          'dt*%s_k2',state_variables);
      case {'rk4','rungekutta','rk'}
        print_k(fid,odes,'_k1',state_variables);                              % write k1 using model.ODEs
        odes_k2=update_odes(odes,'_k1','.5*dt',state_variables,index_lasts);   % F(*,yn+dt*k1/2)
        fprintf(fid,'  t=t+.5*dt;\n');                                          % update t before writing k2
        print_k(fid,odes_k2,'_k2',state_variables);                           % write k2 using odes_k2
        odes_k3=update_odes(odes,'_k2','.5*dt',state_variables,index_lasts);   % F(*,yn+dt*k2/2)
        print_k(fid,odes_k3,'_k3',state_variables);                           % write k3 using odes_k3
        odes_k4=update_odes(odes,'_k3','dt',state_variables,index_lasts);      % F(*,yn+dt*k3)
        fprintf(fid,'  t=t+.5*dt;\n');                                          % update t before writing k4
        print_k(fid,odes_k4,'_k4',state_variables);                           % write k4 using odes_k4
        % write update for state variables
        print_var_update(fid,index_nexts_,index_lasts,...
          '(dt/6)*(%s_k1+2*(%s_k2+%s_k3)+%s_k4)',state_variables);
    end
  end

end
% ########################################
%% SUBFUNCTIONS
% ########################################

function print_k(fid,odes_k,suffix_k,state_variables,nvals_per_var)
  % purpose: write auxiliary calculations (k1-k4) for runge-kutta
  for i=1:length(odes_k)
    fprintf(fid,'  %s%s=%s;\n',state_variables{i},suffix_k,odes_k{i});
  end
end

function odes_out=update_odes(odes,suffix_k,increment,state_variables,index_lasts)
  % purpose: update expressions for axiliary calculations (k1-k4)
  odes_out=odes;
  for i=1:length(odes)
    for j=1:length(odes)
      odes_out{i}=dynasim_strrep(odes_out{i},[state_variables{j} index_lasts{j}],sprintf('(%s%s+%s*%s%s)',state_variables{j},index_lasts{j},increment,state_variables{j},suffix_k),'(',')');    
    end
  end
end

function print_var_update(fid,index_nexts,index_lasts,update_term,state_variables)
  % purpose: write statements to update state variables according to their dynamics
  % example: 
  % update_term='(dt/6)*(%s_k1+2*(%s_k2+%s_k3)+%s_k4)';
  % state_variable='A_v';
  if ~isempty(state_variables)
    fprintf(fid,'  %% ------------------------------------------------------------\n');
    fprintf(fid,'  %% Update state variables:\n');
    fprintf(fid,'  %% ------------------------------------------------------------\n');
  else
    return;
  end
  for i=1:length(state_variables)
    this_update_term=strrep(update_term,'%s',state_variables{i});
    this_update_expression=sprintf('%s%s=%s%s+%s;',state_variables{i},index_nexts{i},state_variables{i},index_lasts{i},this_update_term);
    fprintf(fid,'  %s\n',this_update_expression);
  end
end

function print_var_update_last(fid,index_nexts,state_variables)
  if ~isempty(state_variables)
    fprintf(fid,'  %% ------------------------------------------------------------\n');
    fprintf(fid,'  %% Store state variables:\n');
    fprintf(fid,'  %% ------------------------------------------------------------\n');
  else
    return;
  end
  for i=1:length(state_variables)
    this_update_expression=sprintf('%s%s=%s_last;\n',state_variables{i},index_nexts{i},state_variables{i});
    fprintf(fid,'  %s',this_update_expression);
  end
end

function print_conditional_update(fid,conditionals,index_nexts,state_variables)
  % purpose: write statements to perform conditional actions (that may
  % include updating state variables).
  if ~isempty(conditionals)
    fprintf(fid,'  %% ------------------------------------------------------------\n');
    fprintf(fid,'  %% Conditional actions:\n');
    fprintf(fid,'  %% ------------------------------------------------------------\n');
  else
    return;
  end
  switch index_nexts{1}(1)
    case '('
      action_index='(n,conditional_test)';
    case '_'
      action_index='_last(conditional_test)';    
  end
  for i=1:length(conditionals)
    condition=conditionals(i).condition;
    action=conditionals(i).action;
    elseaction=conditionals(i).else;
    % add indexes to state variables in conditional actions
    for j=1:length(state_variables)
      if strcmp('spike_monitor',conditionals(i).namespace)
        % do nothing if spike_monitor
      else
        condition=dynasim_strrep(condition,state_variables{j},[state_variables{j} index_nexts{j}]);
      end
      action=dynasim_strrep(action,state_variables{j},[state_variables{j} action_index]);
      if ~isempty(elseaction)
        elseaction=dynasim_strrep(elseaction,state_variables{j},[state_variables{j} action_index]);
      end
    end
    % write conditional to solver function
    fprintf(fid,'  conditional_test=(%s);\n',condition);
    action=dynasim_strrep(action,'\(n,:','(n,conditional_test');
    fprintf(fid,'  if any(conditional_test), %s; ',action);
    if ~isempty(elseaction)
      elseaction=dynasim_strrep(elseaction,'(n,:','(n,conditional_test');
      fprintf('else %s; ',elseaction);
    end
    fprintf(fid,'end\n');  
  end
end

function print_monitor_update(fid,monitors,index_nexts,state_variables,index_lasts)
  if ~isempty(monitors) && iscell(index_nexts) % being called from within the integrator loop
    fprintf(fid,'  %% ------------------------------------------------------------\n');
    fprintf(fid,'  %% Update monitors:\n');
    fprintf(fid,'  %% ------------------------------------------------------------\n');
  end
  if nargin<5, index_lasts=index_nexts; end
  % account for inputs from monitor initialization
  if ~iscell(index_nexts), index_nexts={index_nexts}; end
  if ~iscell(index_lasts), index_lasts={index_lasts}; end
  if length(index_lasts)~=length(state_variables)
    index_lasts=repmat(index_lasts,[1 length(state_variables)]); 
  end
  if isequal(index_nexts{1},'(1,:)')
    monitor_index='(1,:)';
  else
    switch index_nexts{1}(1)
      case '('
        monitor_index='(n,:)';
      case '_'
        monitor_index='_last';    
    end  
  end
  % purpose: write statements to update monitors given current state variables
  % note: only run this every time a point is recorded (not necessarily on
  % every step of the integration).
  if isempty(monitors)
    return;
  end
  monitor_name=fieldnames(monitors);
  monitor_expression=struct2cell(monitors);
  % add indexes to state variables in monitors
  for i=1:length(monitor_name)
    for j=1:length(state_variables)
      %monitor_expression{i}=dynasim_strrep(monitor_expression{i},state_variables{j},[state_variables{j} index_lasts{j}]);
      monitor_expression{i}=dynasim_strrep(monitor_expression{i},state_variables{j},[state_variables{j} monitor_index]);
    end
    % write monitors to solver function
    fprintf(fid,'  %s%s=%s;\n',monitor_name{i},monitor_index,monitor_expression{i});
  end
end


%{
PSEUDOCODE:
% --------------------------------------------------------------------------------
% setup (parameters, fixed_variables, functions)
% preallocation and initial conditions
if downsample_factor>1 || disk_flag==1
  % var_last=IC;
  % mon_last=f(IC);
end
if disk_flag==1
  fid=fopen(options.filename,'wt');
  % print var_last
  % print mon_last
else
  % var=zeros(npop,nsamp);
  % mon=zeros(npop,nsamp);
  if downsample_factor==1
    % var(1,:)=IC;
    % mon(1,:)=f(IC);
  else
    % var(1,:)=var_last;
    % mon(1,:)=mon_last;
  end
end  

% numerical integration
% n=1; % storage index
% for k=2:ntime % time index
  if downsample_factor==1 && disk_flag==0 % do not use var_last
    % update_vars;      % var(:,k-1)->var(k,:) or var(k-1)->var(k)
    % conditionals;     % var(k,:)->var(k,:) or var(k)->var(k)
    % update_monitors;  % mon(:,k-1)->mon(k,:) or mon(k-1)->mon(k)
  else % downsample_factor>1 and/or disk_flag==1
    % update_vars;      % var_last->var_last
    % conditionals;     % var_last->var_last
    if mod(t,downsample_factor)==0 % store this time point
      if disk_flag==1 % write to disk
        % write_vars;     % print: var_last
        % write_monitors; % print: mon=f(var_last)
      else            % store in memory
        % update_vars;    % var_last -> var(n,:) or var(n)
        % update_monitors;% f(var_last) -> mon(n,:) or mon(n)
      end
      % n=n+1;
    end
  end
% end
% --------------------------------------------------------------------------------
%}
