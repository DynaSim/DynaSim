function mexfileOutput = PrepareMEX(mfileInput, options)
%PREPAREMEX - take an m-file path and compile it using the Matlab coder.
%
% See also: GetSolveFile

if ~exist(options,'var') || isempty(options)
  options.verbose_flag = 1; % set verbose to 1 by default
end

% make mex name from solve_file
[fpath,fname] = fileparts(mfileInput);
mexfileOutput = fullfile(fpath,[fname '_mex']);

if ~exist(mexfileOutput,'file')
  if options.verbose_flag
    fprintf('Compiling file: %s\n',mfileInput);
    fprintf('            to: %s\n',mexfileOutput);
  end

  compile_start_time=tic;
  
  makeMex(mfileInput); % mex-file solver for solve file

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

end

%% Subfunctions
function makeMex(file)
% Create a MEX configuration object
cfg = coder.config('mex');

% Turn on dynamic memory allocation
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';

% Generate MEX function
eval(['codegen -d codemex -config cfg ',sprintf(file)]);
end
