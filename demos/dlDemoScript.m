%% 3-layer neocortical model of PFC / Predictive task

%% Model parameters, Uncomment one of the following lanes to run an example.

clear;
clc;

Ne = 24;Ni = 6;Nio = 6;noise_rate = 7;
% s0 = NeoCortex(Ne, Ni, Nio, noise_rate);
s1 = dlDemoPING(2, 1, 2, noise_rate);
% s2 = dlDemoPredictivePFC(Ne, Ni, Nio, noise_rate);

%% Create DynaLearn Class (First time)

% m = DynaLearn(s1, 'models/dlTestPING'); % ~ 10min
% m = DynaLearn(s2, 'models/dlDemoPredictivePFC'); % ~ 150min
% m = DynaLearn(s2, 'models/dlDemoPredictivePFC2'); % ~ 180min, 2nd version

m = DynaLearn(s2, 'models/dlTestPredictivePFC', 'raw'); % ~ 42sec
% m.dlSimulate(); % ~ 40sec
m.dlSave(); % < 1sec

%% Load DynaLearn Class (previously saved file is required, default is dlFileBase.mat)

clc;
m = DynaLearn(); % ~ 1sec
% m = m.dlLoad('models/dlDemoPING'); % ~ 10sec
m = m.dlLoad('models/dlDemoPredictivePFC'); % ~ 10sec
% m = m.dlLoad('models/dlDemoPredictivePFC2'); % ~ 10sec
% m = m.dlLoad('models/dlTestPredictivePFC'); % ~ 10sec
% m.dlSimulate(); % ~ 40sec

%% Simulation and general plotting

clc;
Params = containers.Map();
Params('tspan') = [0 1000];
m.dlUpdateParams(Params);

m.dlSimulate(); % (optional) simulate it , ~ seconds runtime
% m.dlPlotAllPotentials('ifr'); % Plot all potential (voltages) as IFR plots ('ifr') or LFP ('lfp').
m.dlPlotAllPotentials('lfp'); % Local field potential

 %% Continue simulation: trialParams example

clc;
trialParams1 = containers.Map();
trialParams2 = containers.Map();
trialParams3 = containers.Map();

trialParams1('tspan') = [0 500];
trialParams2('tspan') = [0 500];
trialParams3('tspan') = [0 500];

g_poisson = 6e-4;dc_poisson = 7e5;

trialParams1('IO_SA1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams1('IO_SA2_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams1('IO_SB1_ctx_iPoisson_DC_poisson') = 0;
trialParams1('IO_SB2_ctx_iPoisson_DC_poisson') = 0;
trialParams1('IO_SC1_ctx_iPoisson_DC_poisson') = 0;
trialParams1('IO_SC2_ctx_iPoisson_DC_poisson') = 0;

trialParams1('IO_SA1_ctx_iPoisson_onset_poisson') = 150;
trialParams1('IO_SA1_ctx_iPoisson_offset_poisson') = 250;
trialParams1('IO_SA2_ctx_iPoisson_onset_poisson') = 250;
trialParams1('IO_SA2_ctx_iPoisson_offset_poisson') = 350;

trialParams1('IO_SB1_ctx_iPoisson_onset_poisson') = 0;
trialParams1('IO_SB1_ctx_iPoisson_offset_poisson') = 0;
trialParams1('IO_SB2_ctx_iPoisson_onset_poisson') = 0;
trialParams1('IO_SB2_ctx_iPoisson_offset_poisson') = 0;

trialParams1('IO_SC1_ctx_iPoisson_onset_poisson') = 0;
trialParams1('IO_SC1_ctx_iPoisson_offset_poisson') = 0;
trialParams1('IO_SC2_ctx_iPoisson_onset_poisson') = 0;
trialParams1('IO_SC2_ctx_iPoisson_offset_poisson') = 0;

trialParams2('IO_SA1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams2('IO_SA2_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams2('IO_SB1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams2('IO_SB2_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams2('IO_SC1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams2('IO_SC2_ctx_iPoisson_DC_poisson') = dc_poisson;

trialParams2('IO_SA1_ctx_iPoisson_onset_poisson') = 250;
trialParams2('IO_SA1_ctx_iPoisson_offset_poisson') = 250;
trialParams2('IO_SA2_ctx_iPoisson_onset_poisson') = 350;
trialParams2('IO_SA2_ctx_iPoisson_offset_poisson') = 350;

trialParams2('IO_SB1_ctx_iPoisson_onset_poisson') = 150;
trialParams2('IO_SB1_ctx_iPoisson_offset_poisson') = 250;
trialParams2('IO_SB2_ctx_iPoisson_onset_poisson') = 250;
trialParams2('IO_SB2_ctx_iPoisson_offset_poisson') = 350;

trialParams2('IO_SC1_ctx_iPoisson_onset_poisson') = 250;
trialParams2('IO_SC1_ctx_iPoisson_offset_poisson') = 250;
trialParams2('IO_SC2_ctx_iPoisson_onset_poisson') = 350;
trialParams2('IO_SC2_ctx_iPoisson_offset_poisson') = 350;

trialParams3('IO_SA1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams3('IO_SA2_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams3('IO_SB1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams3('IO_SB2_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams3('IO_SC1_ctx_iPoisson_DC_poisson') = dc_poisson;
trialParams3('IO_SC2_ctx_iPoisson_DC_poisson') = dc_poisson;

trialParams3('IO_SA1_ctx_iPoisson_onset_poisson') = 250;
trialParams3('IO_SA1_ctx_iPoisson_offset_poisson') = 250;
trialParams3('IO_SA2_ctx_iPoisson_onset_poisson') = 350;
trialParams3('IO_SA2_ctx_iPoisson_offset_poisson') = 350;

trialParams3('IO_SB1_ctx_iPoisson_onset_poisson') = 250;
trialParams3('IO_SB1_ctx_iPoisson_offset_poisson') = 250;
trialParams3('IO_SB2_ctx_iPoisson_onset_poisson') = 350;
trialParams3('IO_SB2_ctx_iPoisson_offset_poisson') = 350;

trialParams3('IO_SC1_ctx_iPoisson_onset_poisson') = 150;
trialParams3('IO_SC1_ctx_iPoisson_offset_poisson') = 250;
trialParams3('IO_SC2_ctx_iPoisson_onset_poisson') = 250;
trialParams3('IO_SC2_ctx_iPoisson_offset_poisson') = 350;

outputParams = [{'DeepE_V', 1:4, [200 400], 'afr'}; {'DeepE_V', 5:8, [200 400], 'afr'}; {'DeepE_V', 9:12, [200 400], 'afr'}; {'DeepE_V', 13:16, [200 400], 'afr'}; {'DeepE_V', 17:20, [200 400], 'afr'}];
targetParams1 = [{'MSE', 1, 25, 0.25}; {'MSE', 2, 12, 0.25}; {'MSE', 3, 12, 0.25}; {'Compare', [1, 2, 3], 0, 0.15}; {'Diff', [2, 3], 0, 0.05}]; % A 
targetParams2 = [{'MSE', 2, 25, 0.25}; {'MSE', 1, 12, 0.25}; {'MSE', 3, 12, 0.25}; {'Compare', [2, 1, 3], 0, 0.15}; {'Diff', [1, 3], 0, 0.05}]; % B
targetParams3 = [{'MSE', 3, 25, 0.25}; {'MSE', 2, 12, 0.25}; {'MSE', 1, 12, 0.25}; {'Compare', [3, 1, 2], 0, 0.15}; {'Diff', [1, 2], 0, 0.05}]; % C

%% Trial: training script 

clc;
dlInputParameters = {trialParams1, trialParams2, trialParams3};
dlTargetParameters = {targetParams1, targetParams2, targetParams3};
dlOutputParameters = outputParams;

dlTrainOptions = containers.Map();
dlTrainOptions('dlEpochs') = 1;
dlTrainOptions('dlBatchs') = 3;
dlTrainOptions('dlLambda') = 1e-5;

dlTrainOptions('dlCheckpoint') = 'true';
dlTrainOptions('dlCheckpointCoefficient') = 1.74; % e.g sqrt(2), sqrt(3), 2, sqrt(5) ... 
dlTrainOptions('dlUpdateMode') = 'batch';
dlTrainOptions('dlLearningRule') = 'NewRule'; % DeltaRule, BioDeltaRule, RWDelta, ...

dlTrainOptions('dlSimulationFlag') = 0; % Manully turning simulation, on or off (on is default and recommended)
dlTrainOptions('dlOutputLogFlag') = 1; % Autosaving trial outputs, on or off (off is default and recommended) % TODO Output/Random/SameValueProblem
dlTrainOptions('dlOfflineOutputGenerator') = 0; % Just for debugging, generates random outputs based on last outputs. 
dlTrainOptions('dlAdaptiveLambda') = 0; % Adaptive lambda parameter; recommended for long simulations.

dlTrainOptions('dlLambdaCap') = 3e-2; % Only if Adaptive lambda is active, recommended to set a upper-bound (UB) or ignore to use default UB (0.01).
% dlTrainOptions('dlMetaLearningRule') = 'true'; % TODOs!

% m.dlResetTraining(); % Reset logs and optimal state error (not the optimal state file)
% m.dlLoadOptimal();  % Load the current optimal state (if exists)
m.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

%% Run a simulation (without training)

for i = 1:1
    m.dlRunSimulation(dlInputParameters{i}, dlOutputParameters);
end

%%
m.dlPlotAllPotentials('lfp');
%%

opts = containers.Map();
% opts("lf") = 50;
% opts("hf") = 100;
% m.dlPlotAllPotentials('avgfft', opts);

opts("lf") = 20;
opts("hf") = 50;
m.dlPlotAllPotentials('avgfft', opts);

%% Errors log plot

clc;
m.dlPlotBatchErrors(3);
% m.dlPlotErrors();

%% End of Demo (7th of May 2022)
