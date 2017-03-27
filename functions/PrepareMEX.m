function mexfileOutput = PrepareMEX(mfileInput, varargin)
%PREPAREMEX - take an m-file path and compile it using the Matlab coder.
%
% See also: GetSolveFile

% Input args
args = varargin;

if isstruct(args{1})
    % User specified an options structure
    keyvals = Options2Keyval(args{1});
    % Convert it to keyvalues and prepend
    keyvals = horzcat(keyvals,args(2:end));
else
    keyvals = args;
end

options=CheckOptions(keyvals,{...
  'verbose_flag',0,{0,1},... % set verbose to 1 by default
  'mexpath',[],[],... % Directory to search for pre-compiled solve files (solve*_mex*)
  },false);

mexpath = options.mexpath;

% make mex name from solve_file
[fpath,fname] = fileparts(mfileInput);
mexfileOutput = fullfile(fpath,[fname '_mex']);

if ~exist(mexfileOutput,'file')
  if options.verbose_flag
    fprintf('Compiling file: %s\n',mfileInput);
    fprintf('            to: %s\n',mexfileOutput);
  end

  compile_start_time=tic;
  
  makeMex(mfileInput, options); % mex-file solver for solve file

  if options.verbose_flag
    fprintf('\tMEX generation complete!\n\tElapsed time: %g seconds.\n',toc(compile_start_time));
    %toc;
  end

  codemex_dir=fullfile(fileparts(mexfileOutput),'codemex');
  if exist(codemex_dir,'dir')
    if options.verbose_flag
      fprintf('\tRemoving temporary codemex directory: %s\n',codemex_dir);
    end
    rmdir(codemex_dir,'s');
  end
else % mex file exists
  if options.verbose_flag
    fprintf('Using previous compiled file: %s\n',mexfileOutput);
  end
end %if

% If mexpath is specified, back up the newly compiled mex files to this
% folder
if ~isempty(mexpath)
    [~,solvefile] = fileparts2(mfileInput);
    [~,mexfile] = fileparts2(mexfileOutput);
    
    if isempty(dir(fullfile(mexpath,[solvefile '.m'])))
        fprintf('Solve file %s does not yet exist in mexpath %s. Copying... \n',solvefile,mexpath);
        if ~exist(mexpath,'dir'); error('Cannot find mexpath! Make sure it exists and is specified as an *absolute* path'); end
        copyfile(mfileInput,mexpath);
    end
    if isempty(dir(fullfile(mexpath,[mexfile '*'])))
        fprintf('Mex file %s does not yet exist in mexpath %s. Copying... \n',mexfile,mexpath);
        if ~exist(mexpath,'dir'); error('Cannot find mexpath! Make sure it exists and is specified as an *absolute* path'); end
        copyfile([mexfileOutput,'*'],mexpath);
    end
end

end

%% Subfunctions
function makeMex(file, options)
% Create a MEX configuration object
cfg = coder.config('mex');

% Turn on dynamic memory allocation
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';

% Generate MEX function
if ~isfield(options, 'codegen_args')
  eval(sprintf('codegen -d codemex -config cfg %s',file));
else % codegen_args specified
  eval(sprintf('codegen -args %s -d codemex -config cfg %s', 'options.codegen_args', file));
end

end
