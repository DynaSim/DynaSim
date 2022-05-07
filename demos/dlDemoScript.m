%% 3-layer neocortical model of PFC / Predictive task

% Manually determined weights version. 
% memoize.m -> suspended due to an error
% dsModel: 

%% Model parameters

clear;
clc;

Ne = 20;Ni = 4;Nio = 10;noise_rate = 13;
% s = NeoCortex(Ne, Ni, Nio, noise_rate);
% s = dlDemoPING(3, 1, 2, noise_rate); % 14 Mins on mex generator
s = dlDemoPredictivePFC(Ne, Ni, Nio, noise_rate);

%% Create DynaLearn Class (First time)

% m = DynaLearn(s, 'models/dlDemoPING'); % ~ 120min
m = DynaLearn(s, 'models/dlDemoPredictivePFC'); % ~ 120min
m.dlSimulate(); % ~ 40sec
m.dlSave(); % < 1sec

%% Load DynaLearn Class (previously saved file is required, default is dlFileBase.mat)

clc;
m = DynaLearn(); % ~ 1sec
% m = m.dlLoad('models/dlDemoPING'); % ~ 10sec
m = m.dlLoad('models/dlDemoPredictivePFC'); % ~ 10sec
m.dlSimulate(); % ~ 40sec

 %% Continue simulation: trialParams example

clc;
g_poisson = 5.7e-4;
dc_poisson = 9e6;

trialParams1 = containers.Map();
trialParams2 = containers.Map();
trialParams3 = containers.Map();

trialParams1('tspan') = [0 500];
trialParams2('tspan') = [0 500];
trialParams3('tspan') = [0 500];

trialParams('SA1_ctx_iPoisson_g_poisson') = g_poisson;
trialParams('SA2_ctx_iPoisson_g_poisson') = g_poisson;
trialParams('SB1_ctx_iPoisson_g_poisson') = g_poisson;
trialParams('SB2_ctx_iPoisson_g_poisson') = g_poisson;
trialParams('SC1_ctx_iPoisson_g_poisson') = g_poisson;
trialParams('SC2_ctx_iPoisson_g_poisson') = g_poisson;

trialParams('SA1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams('SA2_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams('SB1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams('SB2_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams('SC1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams('SC2_ctx_iPoisson_DC_poisson') = dc_poisson;

trialParams('SA1_ctx_iPoisson_onset_poisson') = 150;
trialParams('SA1_ctx_iPoisson_offset_poisson') = 250;
trialParams('SA2_ctx_iPoisson_onset_poisson') = 250;
trialParams('SA2_ctx_iPoisson_offset_poisson') = 350;

trialParams('SB1_ctx_iPoisson_onset_poisson') = 250;
trialParams('SB1_ctx_iPoisson_offset_poisson') = 250;
trialParams('SB2_ctx_iPoisson_onset_poisson') = 350;
trialParams('SB2_ctx_iPoisson_offset_poisson') = 350;

trialParams('SC1_ctx_iPoisson_onset_poisson') = 250;
trialParams('SC1_ctx_iPoisson_offset_poisson') = 250;
trialParams('SC2_ctx_iPoisson_onset_poisson') = 350;
trialParams('SC2_ctx_iPoisson_offset_poisson') = 350;

trialParams1('SA1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams1('SA2_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams1('SB1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams1('SB2_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams1('SC1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams1('SC2_ctx_iPoisson_DC_poisson') = dc_poisson;

trialParams1('SA1_ctx_iPoisson_onset_poisson') = 150;
trialParams1('SA1_ctx_iPoisson_offset_poisson') = 250;
trialParams1('SA2_ctx_iPoisson_onset_poisson') = 250;
trialParams1('SA2_ctx_iPoisson_offset_poisson') = 350;

trialParams1('SB1_ctx_iPoisson_onset_poisson') = 250;
trialParams1('SB1_ctx_iPoisson_offset_poisson') = 250;
trialParams1('SB2_ctx_iPoisson_onset_poisson') = 350;
trialParams1('SB2_ctx_iPoisson_offset_poisson') = 350;

trialParams1('SC1_ctx_iPoisson_onset_poisson') = 250;
trialParams1('SC1_ctx_iPoisson_offset_poisson') = 250;
trialParams1('SC2_ctx_iPoisson_onset_poisson') = 350;
trialParams1('SC2_ctx_iPoisson_offset_poisson') = 350;

trialParams2('SA1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams2('SA2_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams2('SB1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams2('SB2_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams2('SC1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams2('SC2_ctx_iPoisson_DC_poisson') = dc_poisson;

trialParams2('SA1_ctx_iPoisson_onset_poisson') = 250;
trialParams2('SA1_ctx_iPoisson_offset_poisson') = 250;
trialParams2('SA2_ctx_iPoisson_onset_poisson') = 350;
trialParams2('SA2_ctx_iPoisson_offset_poisson') = 350;

trialParams2('SB1_ctx_iPoisson_onset_poisson') = 150;
trialParams2('SB1_ctx_iPoisson_offset_poisson') = 250;
trialParams2('SB2_ctx_iPoisson_onset_poisson') = 250;
trialParams2('SB2_ctx_iPoisson_offset_poisson') = 350;

trialParams2('SC1_ctx_iPoisson_onset_poisson') = 250;
trialParams2('SC1_ctx_iPoisson_offset_poisson') = 250;
trialParams2('SC2_ctx_iPoisson_onset_poisson') = 350;
trialParams2('SC2_ctx_iPoisson_offset_poisson') = 350;

trialParams3('SA1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams3('SA2_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams3('SB1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams3('SB2_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams3('SC1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams3('SC2_ctx_iPoisson_DC_poisson') = dc_poisson;

trialParams3('SA1_ctx_iPoisson_onset_poisson') = 250;
trialParams3('SA1_ctx_iPoisson_offset_poisson') = 250;
trialParams3('SA2_ctx_iPoisson_onset_poisson') = 350;
trialParams3('SA2_ctx_iPoisson_offset_poisson') = 350;

trialParams3('SB1_ctx_iPoisson_onset_poisson') = 250;
trialParams3('SB1_ctx_iPoisson_offset_poisson') = 250;
trialParams3('SB2_ctx_iPoisson_onset_poisson') = 350;
trialParams3('SB2_ctx_iPoisson_offset_poisson') = 350;

trialParams3('SC1_ctx_iPoisson_onset_poisson') = 150;
trialParams3('SC1_ctx_iPoisson_offset_poisson') = 250;
trialParams3('SC2_ctx_iPoisson_onset_poisson') = 250;
trialParams3('SC2_ctx_iPoisson_offset_poisson') = 350;

outputParams = [{'DeepE_V', 1:4, [200 400], 'afr'}; {'DeepE_V', 5:8, [200 400], 'afr'}; {'DeepE_V', 9:12, [200 400], 'afr'}; {'DeepE_V', 13:16, [200 400], 'afr'}; {'DeepE_V', 17:20, [200 400], 'afr'}];
targetParams1 = [{'MSE', 1, 6, 0.25}; {'MSE', 2, 3, 0.25}; {'MSE', 3, 3, 0.25}; {'Compare', [1, 2, 3], 0, 0.15}; {'Diff', [2, 3], 0, 0.05}]; % A 
targetParams2 = [{'MSE', 2, 6, 0.25}; {'MSE', 1, 3, 0.25}; {'MSE', 3, 3, 0.25}; {'Compare', [2, 1, 3], 0, 0.15}; {'Diff', [1, 3], 0, 0.05}]; % B
targetParams3 = [{'MSE', 3, 6, 0.25}; {'MSE', 2, 3, 0.25}; {'MSE', 1, 3, 0.25}; {'Compare', [3, 1, 2], 0, 0.15}; {'Diff', [1, 2], 0, 0.05}]; % C

%% Trial: training script 
% TODO ->>> (similar inputs-outputs problem)

clc;
dlInputParameters = {trialParams1, trialParams2, trialParams3};
dlTargetParameters = {targetParams1, targetParams2, targetParams3};
dlOutputParameters = outputParams;

dlTrainOptions = containers.Map();
dlTrainOptions('dlEpochs') = 2;
dlTrainOptions('dlBatchs') = 3;
dlTrainOptions('dlLambda') = 0.001;

dlTrainOptions('dlCheckpoint') = 'true';
dlTrainOptions('dlCheckpointCoefficient') = 2; % e.g sqrt(2), sqrt(3), 2, sqrt(5) ... 
dlTrainOptions('dlUpdateMode') = 'batch';
dlTrainOptions('dlLearningRule') = 'BioDeltaRule'; % DeltaRule, BioDeltaRule, RWDelta, ...

dlTrainOptions('dlSimulationFlag') = 1; % Manully turning simulation, on or off (on is default and recommended)
dlTrainOptions('dlOutputLogFlag') = 1; % Autosaving trial outputs, on or off (off is default and recommended) % TODO Output/Random/SameValueProblem
dlTrainOptions('dlOfflineOutputGenerator') = 0; % Just for debugging, generates random outputs based on last outputs. 
dlTrainOptions('dlAdaptiveLambda') = 1; % Adaptive lambda parameter; recommended for long simulations.

dlTrainOptions('dlLambdaCap') = 3e-2; % Only if Adaptive lambda is active, recommended to set a upper-bound (UB) or ignore to use default UB (0.01).
% dlTrainOptions('dlMetaLearningRule') = 'TODO'; %%% 
% dlTrainOptions() = '';
% dlTrainOptions() = '';

% m.dlResetTraining(); % Reset logs and optimal state error
m.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

%% Run a simulation (without training)

m.dlRunSimulation(dlInputParameters{1}, dlOutputParameters);
m.dlPlotAllPotentials('lfp');
m.dlRunSimulation(dlInputParameters{2}, dlOutputParameters);
m.dlPlotAllPotentials('lfp');
m.dlRunSimulation(dlInputParameters{3}, dlOutputParameters);
m.dlPlotAllPotentials('lfp');

%% End of Demo