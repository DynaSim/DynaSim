function solve_file_m=CompareSolveFiles(solve_file_m,mexpath)
%COMPARESOLVEFILES - look for an equivalent file in same directory
%
%  - Step 1: compare to other *.m in /solve/
%  - Step 2: if match: remove(solve_file); solve_file=match;
% 
% If mexpath is specified, then will do a comparison to other files
% in the mexpath, as opposed to the current solve folder
%
% See also: GetSolveFile, SimulateModel, CreateBatch

[solvepath,fname,fext]=fileparts(solve_file_m);

if nargin < 2
    mexpath = solvepath;
end

% get list of files in where solve_file is located
D=dir(mexpath);
files={D(~[D.isdir]).name};
files=files(cellfun(@any,regexp(files,'.m$')));
files=setdiff(files,[fname fext]);

% compare solve_file_m to each file
for f=1:length(files)
  [~,diffs] = system(['diff ' solve_file_m ' ' fullfile(mexpath,files{f})]);
  if isempty(diffs)
    %dbstack
    old_solve_file_m=solve_file_m;
    delete(old_solve_file_m);
    
    if ~strcmp(mexpath,solvepath)
        % Copy mex file and solve file from mexpath into solve path
        copyfile(fullfile(mexpath,files{f}),solvepath);
        [~,fname] = fileparts(files{f});
        if ~isempty(dir(fullfile(mexpath,[fname '_mex*'])))
            copyfile(fullfile(mexpath,[fname '_mex*']),solvepath);
        end
    end
    
    solve_file_m=fullfile(solvepath,files{f});
    break;
  end
end
end
