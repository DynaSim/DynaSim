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
  'ic',[],[],...                  % initial conditions (overrides definition in model structure)
  'tspan',[0 100],[],...          % [beg,end] (units must be consistent with dt and equations)
  'dt',.01,[],...                 % time step used for fixed step DynaSim solvers
  'downsample_factor',1,[],...    % downsampling applied after simulation (only every downsample_factor-time point is returned)
  'random_seed','shuffle',[],...        % seed for random number generator (usage: rng(random_seed))
  'solver','euler',{'ode23','ode45','ode1','ode2','ode3','ode4','ode5','ode8',...
    'ode113','ode15s','ode23s','ode23t','ode23tb'},... % DSSim and built-in Matlab solvers
  'matlab_solver_options',[],[],... % options from odeset for use with built-in Matlab solvers
  'reduce_function_calls_flag',1,{0,1},...   % whether to eliminate internal (anonymous) function calls
  'save_parameters_flag',1,{0,1},...
  'solver_type','matlab',{'matlab', 'matlab_no_mex'},... % if compile_flag==1, will decide whether to mex solve_file or odefun_file
  'filename',[],[],...         % name of solver file that integrates model
  'fileID',1,[],...
  'compile_flag',exist('codegen')==6,{0,1},... % whether to prepare script for being compiled using coder instead of interpreting Matlab
  'verbose_flag',1,{0,1},...
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


fprintf(fid,'function %s=solve_ode\n',output_string);

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


%% ODEFUN
if options.compile_flag && strcmp(solver_type, 'matlab_no_mex')
  % save ode function as separate m-file for mex compilation
  
  %open file
  odefun_file = [options.filename 'odefun'];
  odefun_fid=fopen(odefun_file,'wt');
  odefun_outfile=fopen(odefun_fid);
  
  %write to file
  % TODO
  
  %close file
  fclose(odefun_fid);
else
  % make sub function (no shared variables with main function workspace for max performance)
  
  odefun_file = [];
end

if ~strcmp(outfile,'"stdout"')
  fclose(fid);
  % wait for file before continuing to simulation
  while ~exist(outfile,'file')
    pause(.01);
  end
end