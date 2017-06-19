function solve_file_m = compareSolveFiles(solve_file_m,mexPath)
%COMPARESOLVEFILES - look for an equivalent file in same directory
%
%  - Step 1: compare to other *.m in /solve/
%  - Step 2: if match: remove(solve_file); solve_file=match;
% 
% If mexPath is specified, then will do a comparison to other files
% in the mexPath, as opposed to the current solve folder
%
% See also: dsGetSolveFile, dsSimulate, dsCreateBatch

[solvePath,fname,fext]=fileparts(solve_file_m);

if nargin < 2
    mexPath = solvePath;
end

% get list of files in where solve_file is located
D=dir(mexPath);
files={D(~[D.isdir]).name};
files=files(cellfun(@any,regexp(files,'.m$')));
files=setdiff(files,[fname fext]);

% compare solve_file_m to each file
for f=1:length(files)
  [~,diffs] = system(['diff ' solve_file_m ' ' fullfile(mexPath,files{f})]);
  if isempty(diffs)
    %dbstack
    old_solve_file_m=solve_file_m;
    delete(old_solve_file_m);
    
    if ~strcmp(mexPath,solvePath)
        % Copy mex file and solve file from mexPath into solve path
        copyfile(fullfile(mexPath,files{f}),solvePath);
        [~,fname] = fileparts(files{f});
        if ~isempty(dir(fullfile(mexPath,[fname '_mex*'])))
            copyfile(fullfile(mexPath,[fname '_mex*']),solvePath);
        end
    end
    
    solve_file_m=fullfile(solvePath,files{f});
    break;
  end
end
end
