function solve_file_mex = PrepareMEX(solve_file)
%PREPAREMEX - take a solver m-file path and compile it using the Matlab coder.
%
% See also: GetSolveFile

% make mex name from solve_file
[fpath,fname] = fileparts(solve_file);
solve_file_mex=fullfile(fpath,[fname '_mex']);

if ~exist(solve_file_mex,'file')
  if options.verbose_flag
    fprintf('Compiling solver file: %s\n',solve_file_mex);
  end

  compile_start_time=tic;
  
  makeMex(solve_file); % mex-file solver for solve file

  if options.verbose_flag
    fprintf('\tMEX generation complete!\n\tElapsed time: %g seconds.\n',toc(compile_start_time));
    %toc;
  end

  codemex_dir=fullfile(fileparts(solve_file_mex),'codemex');
  if exist(codemex_dir,'dir')
    if options.verbose_flag
      fprintf('\tRemoving temporary codemex directory: %s\n',codemex_dir);
    end
    rmdir(codemex_dir,'s');
  end
else % mex file exists
  if options.verbose_flag
    fprintf('Using previous compiled solver file: %s\n',solve_file_mex);
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