function solve_file_m=CompareSolveFiles(solve_file_m)
% Purpose: look for an equivalent file in same directory
%  - compare to other *.m in /solve/
%  - if match: remove(solve_file); solve_file=match;  
% See also: GetSolveFile, SimulateModel, CreateBatch

[fpath,fname,fext]=fileparts(solve_file_m);
% get list of files in where solve_file is located
D=dir(fpath);
files={D(~[D.isdir]).name};
files=files(cellfun(@any,regexp(files,'.m$')));
files=setdiff(files,[fname fext]);
% compare solve_file_m to each file
for f=1:length(files)
  [~,diffs] = system(['diff ' solve_file_m ' ' fullfile(fpath,files{f})]);
  if isempty(diffs)
    %dbstack
    old_solve_file_m=solve_file_m;
    delete(old_solve_file_m);
    solve_file_m=fullfile(fpath,files{f});
    break;
  end
end
