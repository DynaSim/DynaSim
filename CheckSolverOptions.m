function options=CheckSolverOptions(options)
% purpose: standardize simulation options appended to params.mat
% - use to achieve consistent params.mat whether created by SimulateModel(),
% WriteDynaSimSolver(), or WriteMatlabSolver().

% standardize and set defaults
keyvals = Options2Keyval(options);
options=CheckOptions(keyvals,{...
  'tspan',[0 100],[],...          % [beg,end] (units must be consistent with dt and equations)
  'downsample_factor',1,[],...    % downsampling applied during simulation (only every downsample_factor-time point is stored in memory or written to disk)
  'random_seed','shuffle',[],...        % seed for random number generator (usage: rng(random_seed))
  'solver','rk4',[],... % DynaSim and built-in Matlab solvers
  'disk_flag',0,[],...            % whether to write to disk during simulation instead of storing in memory
  'dt',.01,[],...                 % time step used for fixed step DynaSim solvers
  'datafile','data.csv',[],... % name of data file if disk_flag=1
  'compile_flag',exist('codegen')==6,[],... % whether to prepare script for being compiled using coder instead of interpreting Matlab
  'verbose_flag',1,[],...
  },false);

field_order={'tspan','downsample_factor','random_seed','solver','disk_flag',...
  'dt','datafile','compile_flag','verbose_flag'};

if options.compile_flag==1
  % <-- copied from WriteDynaSimSolver.m -->
  % todo: make seed string (eg, 'shuffle') from param struct work with coder (options.compile_flag=1)
  % (currently raises error: "String input must be constant")
  % workaround: (shuffle here and get numeric seed for MEX-compatible params.mat)
  rng(options.random_seed);
  options.random_seed=getfield(rng,'Seed');
end

% standardize field order
options=orderfields(options,field_order);
