% Below are 2 different ways to simulate on the cluster.  With larger
% simulations (>100s), array mode should be faster.

%% Setup

% Set path to your copy of the DynaSim toolbox
dynasim_path = fullfile('..');

% add DynaSim toolbox to Matlab path
addpath(genpath(dynasim_path)); % comment this out if already in path

% Set where to save outputs
output_directory = 'outputs';

% move to root directory where outputs will be saved
cd(output_directory);

%% Solve with loop mode
study_dir='study_HH_varyI_cluster_loop'; 
eqns='dv/dt=@current+I; {iNa,iK}';
vary={'','I',[0:10:50]};

if exist(study_dir, 'dir')
  rmdir(study_dir, 's')
end

SimulateModel(eqns,'vary',vary, 'study_dir',study_dir,'save_data_flag',1,...
  'cluster_flag',1,'verbose_flag',1,'qsub_mode','loop');

%% Solve with array mode
study_dir='study_HH_varyI_cluster_array'; 
eqns='dv/dt=@current+I; {iNa,iK}';
vary={'','I',[0:10:50]};

if exist(study_dir, 'dir')
  rmdir(study_dir, 's')
end

SimulateModel(eqns,'vary',vary, 'study_dir',study_dir,'save_data_flag',1,...
  'cluster_flag',1,'verbose_flag',1,'qsub_mode','array');