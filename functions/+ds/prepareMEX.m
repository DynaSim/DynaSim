function mexfileOutput = prepareMEX(mfileInput, varargin)
%PREPAREMEX - take an m-file path and compile it using the Matlab coder.
%
% See also: ds.getSolveFile

% Input args
args = varargin;

if isstruct(args{1})
  % User specified an options structure
  keyvals = ds.options2Keyval(args{1});
  % Convert it to keyvalues and prepend
  keyvals = horzcat(keyvals,args(2:end));
else
  keyvals = args;
end

options=ds.checkOptions(keyvals,{...
  'verbose_flag',0,{0,1},... % set verbose to 1 by default
  'mex_dir',[],[],... % Directory to search for pre-compiled solve files (solve*_mex*)
  },false);

mex_dir = options.mex_dir;

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

% If mex_dir is specified, back up the newly compiled mex files to this folder
if ~isempty(mex_dir)
  [~,solvefile] = fileparts2(mfileInput);
  [~,mexfile] = fileparts2(mexfileOutput);
  
  if isempty(dir(fullfile(mex_dir,[solvefile '.m'])))
    if options.verbose_flag
      fprintf('Solve file %s does not yet exist in mex_dir %s. Copying... \n',solvefile,mex_dir);
    end
    if ~exist(mex_dir,'dir'); error('Cannot find %s! Make sure it exists and is specified as an *absolute* path',mex_dir); end
    copyfile(mfileInput,mex_dir);
  end
  
  if isempty(dir(fullfile(mex_dir,[mexfile '*'])))
    if options.verbose_flag
      fprintf('Mex file %s does not yet exist in mex_dir %s. Copying... \n',mexfile,mex_dir);
    end
    if ~exist(mex_dir,'dir'); error('Cannot find %s! Make sure it exists and is specified as an *absolute* path',mex_dir); end
    copyfile([mexfileOutput,'*'],mex_dir);
  end
end

end % main fn

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
