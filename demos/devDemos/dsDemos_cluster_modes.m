% The tests some of the major cluster options

%% Setup

% Add DynaSim to path if it's not already there
if exist('setupDynaSimPath','file')
    setupDynaSimPath;
else
    error('Add the DynaSim folder to the MATLAB path - e.g. run addpath(genpath(DynaSimPath))');
end

% Set where to save outputs
output_directory = dsGetConfig('demos_path');

% move to root directory where outputs will be saved
mkdirSilent(output_directory);
cd(output_directory);

%% options
qsub_modes = {'loop', 'array'};
% compile_flag = [0,1];
% one_solve_file_flag = [0,1];
% parallel_flag = [0,1];

N = 2;
v = (1:2^(N^2))-1;
opts = dec2bin(v)' - '0';
nOpt = size(opts,2);
for iOpt = 1:nOpt
  thisOpt = opts(:, iOpt);
  
  % get opt values
  qsub_mode = qsub_modes{thisOpt(1)+1};
  compile_flag = thisOpt(2);
  parallel_flag = thisOpt(3);
  one_solve_file_flag = thisOpt(4);
  
  % make dir name
  study_dir = fullfile(output_directory, sprintf('study_HH_varyI_cluster_%s',qsub_mode));
  if compile_flag
    study_dir = [study_dir '_comp'];
  end
  if parallel_flag
    study_dir = [study_dir '_par'];
  end
  if one_solve_file_flag
    study_dir = [study_dir '_oneFile'];
  end

  dsSimulate(eqns,'vary',vary, 'study_dir',study_dir,'save_data_flag',1,...
  'cluster_flag',1,'verbose_flag',1,'qsub_mode',qsub_mode, 'compile_flag',compile_flag,...
  'one_solve_file_flag',one_solve_file_flag, 'parallel_flag',parallel_flag);
end