% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DynaSim Classification Demos
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Download DynaSim toolbox from https://github.com/dynasim/dynasim.
  or, download using git: git clone https://github.com/dynasim/dynasim.git
For further documentation, see tutorial.m in the demos directory
Sign up for user mailing list at: https://groups.google.com/forum/#!forum/dynasim-users.

Tip: In Matlab, you can obtain more information associated with any function "FUNCTION_NAME"
by entering "help FUNCTION_NAME" in the command window. Use the "See also" list
at the end of the help section to browse through related help documentation.
%}

% Get ready...

% Set path to your copy of the DynaSim toolbox
dynasim_path = '..';

% add DynaSim toolbox to Matlab path
addpath(genpath(dynasim_path)); % comment this out if already in path

% Set where to save outputs
study_dir = fullfile(pwd, 'outputs', 'demo_sPING_classify');

% Here we go!

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILDING LARGE MODELS WITH MULTIPLE POPULATIONS AND CONNECTIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sparse Pyramidal-Interneuron-Network-Gamma (sPING)

% define equations of cell model (same for E and I populations)
eqns={
  'dV/dt=Iapp+@current+noise*randn(1,N_pop)';
  };
% Tip: monitor all functions of a mechanism using: monitor MECHANISM.functions

% create DynaSim specification structure
s=[];
% neural pops
s.populations(1).name='E';
s.populations(1).size=80;
s.populations(1).equations=eqns;
s.populations(1).mechanism_list={'iNa','iK'};
s.populations(1).parameters={'Iapp',0,'gNa',120,'gK',36,'noise',40};
s.populations(2).name='I';
s.populations(2).size=20;
s.populations(2).equations=eqns;
s.populations(2).mechanism_list={'iNa','iK'};
s.populations(2).parameters={'Iapp',0,'gNa',120,'gK',36,'noise',40};

% synapses
s.connections(1).direction='I->E';
s.connections(1).mechanism_list={'iGABAaScaled'};
s.connections(1).parameters={'tauD',10, 'gSYN',1, 'prob_cxn',1};
s.connections(2).direction='E->I';
s.connections(2).mechanism_list={'iAMPAScaled'};
s.connections(2).parameters={'tauD',2, 'gSYN',1, 'prob_cxn',1};

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVING SIMULATED DATA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How to: set 'save_data_flag' to 1 and (optionally) 'study_dir' to /path/to/outputs

%% Save data from a set of simulations
time_end = 500;
dt = 0.01;

% Specify what to vary
% Tip: use 'vary' Syntax to to systematically vary a parameter
vary={
  'E',     'Iapp',  [0:0.5:1];
  'I',     'Iapp',  [0:0.5:1];
%   '(E,I)', 'Noise', [0:10:30];
%   'E',     'prob_cxn', [1, 0.5, .1];
%   'I',     'prob_cxn', [1, 0.5, .1];
  'E',     'gSYN', [1, 0.5];
  'I',     'gSYN', [1, 0.5];
  }; % vary the amplitude of tonic input to E-cells

% if exist(study_dir,'dir')
%   rmdir(study_dir,'s')
% end

data=SimulateModel(s,'save_data_flag',0, 'save_results_flag',1, 'overwrite_flag',0, 'study_dir',study_dir, 'compile_flag',1,...
  'vary',vary, 'solver','euler', 'dt',0.01, 'verbose_flag',1, 'tspan', [0 time_end], 'downsample_factor',1/dt,...
  'analysis_functions',{@classifyEI},...
  'plot_functions',{@PlotData,@PlotData},...
  'plot_options',{
    {'varied_filename_flag', 1, 'format', 'jpg', 'visible', 'off'},...
    {'plot_type','rastergram','varied_filename_flag', 1, 'format', 'jpg', 'visible', 'off'},...
  });

%% load and plot the saved data
study_dir = fullfile(pwd, 'outputs', 'demo_sPING_classify');
% if ~exist('data','var')
%   data=ImportData('demo_sPING_classify');
% end
% PlotData(data);
% PlotData(data,'plot_type','rastergram');

% xAxisVaryParamInd = 2;
% yAxisVaryParamInd = 1;
% PlotClass(study_dir, xAxisVaryParamInd, yAxisVaryParamInd);

% gvRunDS(study_dir, struct('overwrite',1))
gvRunDS(study_dir, struct('overwrite',0))