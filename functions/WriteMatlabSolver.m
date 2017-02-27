function data=WriteMatlabSolver(model,varargin)
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
% See also: SimulateModel, DynaSim2Odefun

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



