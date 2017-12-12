function solve_file_m = dsCompareSolveFiles(solve_file_m,mexPath,verbose_flag)
%COMPARESOLVEFILES - look for an equivalent file in same directory
%
%  - Step 1: compare to other *.m in /solve/
%  - Step 2: if match: remove(solve_file); solve_file=match;
% 
% If mexPath is specified, then will do a comparison to other files
% in the mexPath, as opposed to the current solve folder
%
% See also: dsGetSolveFile, dsSimulate, dsCreateBatch
%
% % Example code for testing mex_dir options:
% eqns={
%   's=10; r=27; b=2.666'
%   'dx/dt=s*(y-x)'
%   'dy/dt=r*x-y-x*z'
%   'dz/dt=-b*z+x*y'
% };
% data=dsSimulate(eqns, 'tspan',[0 100], 'ic',[1 2 .5],'verbose',1, 'solver','rk4', 'study_dir','demo_lorenz','compile_flag',1,'mex_dir_flag',0,'mex_dir',[]);
% data=dsSimulate(eqns, 'tspan',[0 100], 'ic',[1 2 .5],'verbose',1, 'solver','rk4', 'study_dir','demo_lorenz','compile_flag',1,'mex_dir_flag',1,'mex_dir',[]);
% data=dsSimulate(eqns, 'tspan',[0 100], 'ic',[1 2 .5],'verbose',1, 'solver','rk4', 'study_dir','demo_lorenz','compile_flag',1,'mex_dir_flag',1,'mex_dir','mexes_temp');
% 
% Author: Jason Sherfey, PhD <jssherfey@gmail.com>
% Copyright (C) 2016 Jason Sherfey, Boston University, USA

[solvePath,fname,fext]=fileparts(solve_file_m);

if nargin < 2
    mexPath = solvePath;
end

if nargin < 3
    verbose_flag = false;
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
            
            if verbose_flag
                fprintf('Copying mex file from %s to %s \n', fullfile(mexPath,[fname '_mex*']), solvePath);
            end
        end
    end
    
    solve_file_m=fullfile(solvePath,files{f});
    break;
  end
end
end
