function [outfile, odefun_file]=WriteMatlabSolver(model,varargin)
%WRITEMATLABSOLVER - write m-file that numerically inteegrates the model
%
% Usage:
%   outfile=WriteMatlabSolver(model,varargin)
%
% Inputs:
%   - model: DynaSim model structure (see GenerateModel)
%   - options:
%     'tspan'         : units must be consistent with dt and equations
%                       {[beg,end]} (default: [0 100])
%     'ic'            : initial conditions; this overrides definition in model structure
%     'solver'        : DynaSim and built-in Matlab solvers
%                       {'euler','rk2','rk4','modifiedeuler',
%                       'ode23','ode45','ode15s','ode23s'}
%     'solver_options': options from odeset for use with built-in Matlab solvers
%     'dt'            :  time step used for fixed step DSSim solvers (default: 0.01)
%     'modifications' : DynaSim modifications structure
%     'reduce_function_calls_flag': whether to eliminate internal function
%                                   calls {0 or 1} (default: 1)
%     'coder_flag'    : whether to compile using coder instead of interpreting
%                       Matlab (default: exist('codegen')==6 TODO is this correct?
%                       what does this mean?)
%     'disk_flag'     : whether to write to disk during simulation instead of
%                       storing in memory {0 or 1} (default: 0)
%     'downsample_factor': downsampling applied during simulation. Only every
%                          downsample_factor-time point is stored in memory or
%                          written to disk (default: 1)
%     'random_seed'   : seed for random number generator (usage:
%                       rng(random_seed)) (default: now)
%
% Outputs:
%   - outfile (solve_ode.m)
%
% Dependencies: CheckOptions, CheckModel
%
% See also: SimulateModel, dynasim2odefun

% Check inputs
options=CheckOptions(varargin,{...
  'tspan',[0 100],[],...          % [beg,end] (units must be consistent with dt and equations)
  'ic',[],[],...                  % initial conditions (overrides definition in model structure)
  'downsample_factor',1,[],...    % downsampling applied after simulation (only every downsample_factor-time point is returned)
  'random_seed','shuffle',[],...        % seed for random number generator (usage: rng(random_seed))
  'solver','euler',{'ode23','ode45'},... % DSSim and built-in Matlab solvers
  'matlab_solver_options',[],[],... % options from odeset for use with built-in Matlab solvers
%   'compile_flag',0,{0,1},... % whether to compile using coder instead of interpreting Matlab  
  'compile_flag',exist('codegen')==6,{0,1},... % whether to prepare script for being compiled using coder instead of interpreting Matlab
  'solver_type','matlab',{'matlab', 'matlab_no_mex'},... % if compile_flag==1, will decide whether to mex solve_file or odefun_file
  'filename',[],[],...         % name of solver file that integrates model
  },false);

% Check inputs
model=CheckModel(model); 

% convert matlab solver options from key/value to struct using odeset if necessary
if iscell(options.matlab_solver_options) && ~isempty(options.matlab_solver_options)
  options.matlab_solver_options=odeset(options.matlab_solver_options{:});
end

%% 1.0 Get ode_fun

% create function that calls feval(@solver,...) and has subfunction
% defining odefun (including optional conditionals)...

[odefun,IC,elem_names]=dynasim2odefun(PropagateParameters(PropagateFunctions(model)), 'odefun_output','func_body');
keyboard
%% 2.0 write m-file solver
% 2.1 create mfile
if ~isempty(options.filename)
  if options.verbose_flag
    fprintf('Creating solver file: %s\n',options.filename);
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

%%
% evaluate fixed_variables
% ...
% solve
if ~isempty(options.solver_options)
  [time,data]=feval(parms.solver,odefun,options.tspan,options.ic,options.solver_options);
else
  [time,data]=feval(parms.solver,odefun,options.tspan,options.ic);
end

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

%% save ode function as separate for mex compilation if needed
if options.compile_flag && strcmp(solver_type, 'matlab_no_mex')
  %save odefun_file
end