%% Simulation of predictive task implemented on DynaSim/DynaLearn via Reinforcement learning 

% Based on [Bastos2020:Layer and rhythm specificity for predictive routing, A.M Bastos et. al, 2020]
% Requirement: DynaSim/dev latest version
% Minimum MatLab version: R2019a
% 5th of July, 2022
% @HNXJ

%% Model parameters

clear;clc;
Ne = 24;Ni = 6;Nin = 6;noise_rate = 5;
s3 = dlModelPredictivePFC(Ne, Ni, Nin, noise_rate); % Predictive PFC model with specific parameters

%% Create DynaLearn Class (Only first time, if file does not exist already)

% m = DynaLearn(s3, 'models/dlModelPredictivePFC2'); % ~180 min, MEXGEN
m.dlSave(); % < 1sec

%% Load DynaLearn Class

m = DynaLearn(); % ~ 1sec
m = m.dlLoad('models/dlModelPredictivePFC2'); % ~ 10sec
% m.dlSimulate(); % ~ 40sec

 %% Continue simulation: trialParams example

[trialParams1, trialParams2, trialParams3] = dlDemoThreePatternNew();

outputParams = [{'DeepE_V', 1:4, [250 450], 'afr'}; {'DeepE_V', 5:8, [250 450], 'afr'}; {'DeepE_V', 9:12, [250 450], 'afr'}];
targetParams1 = [{'MSE', 1, 21, 0.25}; {'MSE', 2, 14, 0.25}; {'MSE', 3, 14, 0.25}; {'Compare', [1, 2, 3], 0, 0.15}; {'Diff', [2, 3], 0, 0.05}]; % A 
targetParams2 = [{'MSE', 2, 21, 0.25}; {'MSE', 1, 14, 0.25}; {'MSE', 3, 14, 0.25}; {'Compare', [2, 1, 3], 0, 0.15}; {'Diff', [1, 3], 0, 0.05}]; % B
targetParams3 = [{'MSE', 3, 21, 0.25}; {'MSE', 2, 14, 0.25}; {'MSE', 1, 14, 0.25}; {'Compare', [3, 1, 2], 0, 0.15}; {'Diff', [1, 2], 0, 0.05}]; % C

%% Trial: training script preparation, 50-block and 50-trial

dlInputParameters = {trialParams1, trialParams2, trialParams3};
dlTargetParameters = {targetParams1, targetParams2, targetParams3};
dlOutputParameters = outputParams;

[B1, B2, B3, T1, T2, T3, TrB, TrT] = dlTrialBlockGenerator(dlInputParameters, dlTargetParameters, 50, 50);

dlTrainOptions = containers.Map();
dlTrainOptions('dlEpochs') = 1;
dlTrainOptions('dlBatchs') = 1;
dlTrainOptions('dlLambda') = 1e-5;
    
dlTrainOptions('dlCheckpoint') = 'true';
dlTrainOptions('dlCheckpointCoefficient') = 1.94; % e.g sqrt(2), sqrt(3), 2, sqrt(5) ... 
dlTrainOptions('dlUpdateMode') = 'trial';
dlTrainOptions('dlLearningRule') = 'BioDeltaRule'; % DeltaRule, BioDeltaRule, RWDelta, ...

dlTrainOptions('dlSimulationFlag') = 1; % Manully turning simulation, on or off (on is default and recommended)
dlTrainOptions('dlOutputLogFlag') = 1; % Autosaving trial outputs, on or off (off is default and recommended) % TODO Output/Random/SameValueProblem
dlTrainOptions('dlOfflineOutputGenerator') = 0; % Just for debugging, generates random outputs based on last outputs. 
dlTrainOptions('dlAdaptiveLambda') = 0; % Adaptive lambda parameter; recommended for long simulations.

dlTrainOptions('dlLambdaCap') = 3e-2; % Only if Adaptive lambda is active, recommended to set a upper-bound (UB) or ignore to use default UB (0.01).
% dlTrainOptions('dlMetaLearningRule') = 'true'; % TODOs!

%% Train test

% m.dlResetTraining(); % Reset logs and optimal state error (not the optimal state file)
% m.dlLoadOptimal();  % Load the current optimal state (if exists)
m.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);
% Raw simulation on (20*3 = 60) iterations

%% Block phase

m.dlTrain(B1, dlOutputParameters, T1, dlTrainOptions);
m.dlTrain(TrB1, dlOutputParameters, TrT1, dlTrainOptions);
m.dlTrain(B2, dlOutputParameters, T2, dlTrainOptions);
m.dlTrain(TrB1, dlOutputParameters, TrT1, dlTrainOptions);
m.dlTrain(B3, dlOutputParameters, T3, dlTrainOptions);
m.dlTrain(TrB1, dlOutputParameters, TrT1, dlTrainOptions);

%% Errors log plot

% clc;
m.dlPlotBatchErrors(3);

%% Plot Local-field potentials

clc;
m.dlPlotAllPotentials('lfp');

%% Run a simulation (without training)

for i = 1:3
    m.dlRunSimulation(dlInputParameters{i}, dlOutputParameters);
%     m.dlPlotAllPotentials('lfp');
end

%%

opts = containers.Map();
% opts("lf") = 50;
% opts("hf") = 100;
% m.dlPlotAllPotentials('avgfft', opts);

opts("lf") = 20;
opts("hf") = 50;
m.dlPlotAllPotentials('avgfft', opts);

%% End of Demo (14th of June 2022)
% Appendix

clc;
fn1 = fieldnames(a1);
vl1 = struct2cell(a1);

fn2 = fieldnames(a2);
vl2 = struct2cell(a2);

ncons = find(contains(fn1, '_netcon'));
for i = ncons'
    x = size(vl1{i, 1});
    y = size(vl2{i, 1});
    fprintf("%d %d %s %s\n", x(1), x(2), fn1{i, 1}, fn2{i, 1});
end
