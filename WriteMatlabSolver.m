function data=WriteMatlabSolver(model,varargin)
%% outfile=WriteMatlabSolver(model,varargin)
% Purpose: write m-file that numerically inteegrates the model
% inputs:
%   model: DSSim model structure (see GenerateModel)
%   options:
%     'tspan',[0 100],[],...          % [beg,end] (units must be consistent with dt and equations)
%     'ic',[],[],...                  % initial conditions (overrides definition in model structure)
%     'solver','euler',{'euler','rk2','rk4','modifiedeuler','ode23','ode45','ode15s','ode23s'},... % DSSim and built-in Matlab solvers
%     'solver_options',[],[],...      % options from odeset for use with built-in Matlab solvers
%     'dt',.01,[],...                 % time step used for fixed step DSSim solvers
%     'modifications',[],[],...       % *DSSim modifications structure
%     'reduce_function_calls_flag',1,[],...   % whether to eliminate internal function calls
%     'coder_flag',exist('codegen')==6,[],... % whether to compile using coder instead of interpreting Matlab
%     'disk_flag',0,[],...            % whether to write to disk during simulation instead of storing in memory
%     'downsample_factor',1,[],...    % downsampling applied during simulation (only every downsample_factor-time point is stored in memory or written to disk)
%     'random_seed',now,[],...        % seed for random number generator (usage: rng(random_seed))
% outputs:
%   outfile (solve_ode.m)
% 
% see also: SimulateModel, DSSimToOdefun
% dependencies: CheckOptions, CheckModel

% Check inputs
options=CheckOptions(varargin,{...
  'tspan',[0 100],[],...          % [beg,end] (units must be consistent with dt and equations)
  'ic',[],[],...                  % initial conditions (overrides definition in model structure)
  'downsample_factor',1,[],...    % downsampling applied after simulation (only every downsample_factor-time point is returned)
  'random_seed','shuffle',[],...        % seed for random number generator (usage: rng(random_seed))
  'solver','euler',{'ode23','ode45'},... % DSSim and built-in Matlab solvers
  'matlab_solver_options',[],[],... % options from odeset for use with built-in Matlab solvers
  },false);

% Check inputs
model=CheckModel(model); 

% convert matlab solver options from key/value to struct using odeset if necessary
if iscell(options.matlab_solver_options) && ~isempty(options.matlab_solver_options)
  options.matlab_solver_options=odeset(options.matlab_solver_options{:});
end

%% 1.0 ...

odefun=DSSimToOdefun(model);

% create function that calls feval(@solver,...) and has subfunction
% defining odefun (including optional conditionals)...

% evaluate fixed_variables
% ...
% solve
if ~isempty(options.solver_options)
  [time,data]=feval(parms.solver,odefun,options.tspan,options.ic,options.solver_options);
else
  [time,data]=feval(parms.solver,odefun,options.tspan,options.ic);
end



