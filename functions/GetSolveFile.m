function solve_file = GetSolveFile(model,studyinfo,varargin)
%% solve_file = GetSolveFile(model,studyinfo,options)
% Purpose: helper function that creates or retrieves the desired solver file.
% Inputs:
%   model - DynaSim model structure (see GenerateModel)
%   studyinfo (optional) - DynaSim studyinfo structure (see CheckStudyinfo)
%   options (optional) - cell array of key/value pairs or Matlab structure with options
%     'solver'      : solver for numerical integration (see GetSolveFile)
%                     {'euler','rk2','rk4'} (default: 'rk4')
%     'disk_flag'     : whether to write to disk during simulation instead of storing in memory {0 or 1} (default: 0)
%     'study_dir'     : relative or absolute path to output directory (default: current directory)
%     'verbose_flag'  : whether to display informative messages/logs (default: 0)
% 
% Output:
%   solver_file - full file name of file solving the system in model
% 
% See also: WriteDynaSimSolver, CompareSolveFiles, PrepareMEX, 
%           SimulateModel, CreateBatch

% Check inputs
opts=[];
if nargin<1
  error('first argument must be a DynaSim model structure.');
end
if nargin<2
  studyinfo=[];
end
if nargin<3
  varargin={};
elseif isstruct(varargin{1}) % user provided an options structure
  opts=varargin{1};
  fields=fieldnames(opts);
  varargin={};
end
options=CheckOptions(varargin,{...
  'solver','rk4',{'euler','rk1','rk2','rk4','modified_euler','rungekutta','rk','ode23','ode45'},... % DynaSim and built-in Matlab solvers
  'matlab_solver_options',[],[],... % options from odeset for use with built-in Matlab solvers
  'disk_flag',0,{0,1},...            % whether to write to disk during simulation instead of storing in memory
  'solve_file',[],[],... % m- or mex-file solving the system
  'study_dir',[],[],... % study directory
  'verbose_flag',0,{0,1},...
  'parallel_flag',0,{0,1},...     % whether to run simulations in parallel (using parfor)
  'compile_flag',0,{0,1},... % exist('codegen')==6, whether to compile using coder instead of interpreting Matlab  
  },false);
if ~isempty(opts)
  % combine default options and user-supplied options w/ the latter
  % overriding the former
  warning('off','catstruct:DuplicatesFound');
  options=catstruct(options,opts);
  options=orderfields(options,fields);
end

if options.verbose_flag
  fprintf('PREPARING SOLVER:\n');
end

% check solver options
switch options.solver
  case {'euler','rk1','rk2','modified_euler','rk4','rungekutta','rk'}
    solver_type = 'dynasim';
    if ~isempty(options.matlab_solver_options)
      warning('matlab_solver_options are not used by DynaSim solvers. instead try ''ode23'' or ''ode45'' to use those options.');
    end
  case {'ode23','ode45'} % note: only ode23 and ode45 are supported by codegen (see: http://www.mathworks.com/help/coder/ug/functions-supported-for-code-generation--alphabetical-list.html)
    solver_type = 'matlab';
    if options.disk_flag==1
      warning('using disk for real-time storage instead of memory is only available for DynaSim solvers. try using ''euler'',''rk2'', or ''rk4'' for real-time disk usage.');
    end
  otherwise
    error('unrecognized solver type');
end

if ~isempty(options.solve_file)
  % use user-provided solve_file
  solve_file=options.solve_file;
  % note: options.solve_file is used by cluster sim jobs (see CreateBatch())
elseif isfield(studyinfo,'solve_file')
  % use study-associated solve_file
  solve_file=studyinfo.solve_file;
else
  % set default solve_file name
  solve_file=['solve_ode_' datestr(now,'yyyymmddHHMMSS_FFF') '.m'];
end
[fpath,fname,fext]=fileparts(solve_file);
if isempty(fpath)
  % add path to solve_file name
  if ~isempty(options.sim_id)
    solve_file=fullfile(options.study_dir,'solve',['sim' num2str(options.sim_id)],[fname fext]);
  else
    solve_file=fullfile(options.study_dir,'solve',[fname fext]);
  end
  % convert relative path to absolute path
  if ~strcmp('/',solve_file(1)) && ~strcmp('\',solve_file(1))
    solve_file=fullfile(pwd,solve_file);
  end
end
[fpath,fname,fext]=fileparts(solve_file);
% check that solve file name is less than max function name allwoed by matlab
if length(fname)>(63-4) % subtract 4 to allow suffix '_mex'
  fname=fname(1:(63-4));
  solve_file=fullfile(fpath,[fname fext]);
end
% create directory for solve_file if it doesn't exist
if ~isdir(fpath)
  if options.verbose_flag
    fprintf('creating solver directory %s\n',fpath);
  end
  mkdir(fpath);
end
cwd=pwd;
if ~strcmp(cwd,fpath)
  if options.verbose_flag
    fprintf('changing directory to %s\n',fpath);
  end
  cd(fpath);
end
% create solve_file if it doesn't exist
if ~exist(solve_file,'file')
  keyvals = Options2Keyval(options);
  switch solver_type
    case 'dynasim'  % write DynaSim solver function (solve_ode.m)
      solve_file_m = WriteDynaSimSolver(model,keyvals{:},'filename',solve_file); % create DynaSim solver m-file
    case 'matlab' % prepare model function handle etc for built-in solver (@odefun)
      solve_file_m = WriteMatlabSolver(model,keyvals{:},'filename',solve_file); % create Matlab solver m-file
                % design: WriteMatlabSolver should be very similar to
                % WriteDynaSimSolver except have a subfunction with an odefun
                % format variation and main function that calls odeset and
                % feval. 
                % DynaSimToOdefun(): a function called outside of
                % SimulateModel. it should evaluate fixed_variables and
                % return @odefun with all substitutions. SimulateModel
                % should be able to handle: SimulateModel(@odefun,'tspan',tspan,'ic',ic)
  end
  solve_file=CompareSolveFiles(solve_file_m);
else
  if options.verbose_flag
    fprintf('using previous solver file: %s\n',solve_file);
  end      
end
% create MEX file if desired and doesn't exist
[fpath,fname,fext]=fileparts(solve_file);
solve_file_mex=fullfile(fpath,[fname '_mex']);
if options.compile_flag % compile solver function        
  if ~exist(solve_file_mex,'file')
    if options.verbose_flag
      fprintf('compiling solver file: %s\n',solve_file_mex);
    end
    compile_start_time=tic;
    PrepareMEX(solve_file); % mex-file solver
    if options.verbose_flag
      fprintf('\tMEX generation complete!\n\tElapsed time: %g seconds.\n',toc(compile_start_time));
      %toc;
    end
    codemex_dir=fullfile(fileparts(solve_file_mex),'codemex');
    if exist(codemex_dir,'dir')
      if options.verbose_flag
        fprintf('\tremoving temporary codemex directory: %s\n',codemex_dir);
      end
      rmdir(codemex_dir,'s');
    end
  else
    if options.verbose_flag
      fprintf('using previous compiled solver file: %s\n',solve_file_mex);
    end
  end
  solve_file=solve_file_mex;
end
if ~strcmp(cwd,fpath)
  if options.verbose_flag
    fprintf('changing directory back to %s\n',cwd);
  end
  cd(cwd);
end