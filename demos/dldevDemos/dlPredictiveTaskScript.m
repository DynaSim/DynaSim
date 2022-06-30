%% Simulation of predictive task implemented on DynaSim/DynaLearn via Reinforcement learning 

% Based on [Bastos2020:Layer and rhythm specificity for predictive routing, A.M Bastos et. al, 2020]
% Requirement: DynaSim/dev latest version
% Minimum MatLab version: R2019a
% 5th of July, 2022
% @HNXJ

%% Model parameters

clear;clc;
Ne = 24;Ni = 6;Nin = 7;NoiseRate = 5;
s3 = dlModelPredictivePFC(Ne, Ni, Nin, NoiseRate); % Predictive PFC model with specific parameters

%% Create DynaLearn Class (Only first time, if file does not exist already)

% m = DynaLearn(s3, 'models/dlModelPredictivePFC3'); % ~70 min, MEXGEN
% m.dlSave(); % < 1sec

%% Load DynaLearn Class

m = DynaLearn(); % ~ 1sec
m = m.dlLoad('models/dlModelPredictivePFC3'); % ~ 10sec
% m.dlSimulate(); % ~ 40sec

 %% Continue simulation: trialParams example

[trialParams1, trialParams2, trialParams3] = dlDemoThreePattern();

outputParams = [{'DeepE_V', 1:4, [300 500], 'afr'}; {'DeepE_V', 5:8, [300 500], 'afr'}; {'DeepE_V', 9:12, [300 500], 'afr'}];
targetParams1 = [{'MSE', 1, 24, 0.2}; {'MSE', 2, 15, 0.2}; {'MSE', 3, 15, 0.2}; {'Compare', [1, 2], 0, 0.3}; {'Compare', [1, 3], 0, 0.3}; {'Diff', [2, 3], 0, 0.04}]; % A 
targetParams2 = [{'MSE', 2, 24, 0.2}; {'MSE', 1, 15, 0.2}; {'MSE', 3, 15, 0.2}; {'Compare', [2, 1], 0, 0.3}; {'Compare', [2, 3], 0, 0.3}; {'Diff', [1, 3], 0, 0.04}]; % B
targetParams3 = [{'MSE', 3, 24, 0.2}; {'MSE', 2, 15, 0.2}; {'MSE', 1, 15, 0.2}; {'Compare', [3, 1], 0, 0.3}; {'Compare', [3, 2], 0, 0.3}; {'Diff', [1, 2], 0, 0.04}]; % C

%% Trial: training script preparation, 50-block and 50-trial

dlInputParameters = {trialParams1, trialParams2, trialParams3};
dlTargetParameters = {targetParams1, targetParams2, targetParams3};
dlOutputParameters = outputParams;

TBdata = dlTrialBlockGenerator(dlInputParameters, dlTargetParameters, 50, 50);

dlTrainOptions = containers.Map();
dlTrainOptions('dlEpochs') = 500;
dlTrainOptions('dlBatchs') = 3;
dlTrainOptions('dlLambda') = 1e-5;
    
dlTrainOptions('dlCheckpoint') = 'true';
dlTrainOptions('dlCheckpointCoefficient') = 1.74; % A.K.A exploration rate 
dlTrainOptions('dlUpdateMode') = 'batch';
dlTrainOptions('dlLearningRule') = 'BioDeltaRule'; % Delta rule with a basic change based on biophysical properties 

dlTrainOptions('dlSimulationFlag') = 1;
dlTrainOptions('dlOutputLogFlag') = 1;
dlTrainOptions('dlOfflineOutputGenerator') = 0;
dlTrainOptions('dlAdaptiveLambda') = 0; % Adaptive lambda parameter; recommended for long simulations.

dlTrainOptions('dlLambdaCap') = 3e-2; % Only if Adaptive lambda is active, recommended to set a upper-bound (UB) or ignore to use default UB (0.01).

%% Train test
% Initial training on the model to reach a plausible local minimia like
% the task in the paper the model should also learn the basics of the task.
% We shortly train the model by cues to put it close to a local minimia.

dlTrainOptions('dlCheckpointCoefficient') = 1.74;
dlTrainOptions('dlLambda') = 6e-6;
dlTrainOptions('dlEpochs') = 1000;
dlTrainOptions('dlBatchs') = 3;

m.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

%% Block-trial phase

clc;
dlTrainOptions('dlLambda') = 3e-6;
dlTrainOptions('dlUpdateMode') = 'trial';
dlTrainOptions('dlEpochs') = 1;
dlTrainOptions('dlBatchs') = 50;

dlTrainOptions('dlCheckpointCoefficient') = 4;
m.dlTrain(TBdata.B1, dlOutputParameters, TBdata.T1, dlTrainOptions);
dlTrainOptions('dlCheckpointCoefficient') = 7;
m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);

dlTrainOptions('dlCheckpointCoefficient') = 4; 
m.dlTrain(TBdata.B2, dlOutputParameters, TBdata.T2, dlTrainOptions);
dlTrainOptions('dlCheckpointCoefficient') = 7; 
m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);

dlTrainOptions('dlCheckpointCoefficient') = 4;
m.dlTrain(TBdata.B3, dlOutputParameters, TBdata.T3, dlTrainOptions);
dlTrainOptions('dlCheckpointCoefficient') = 7;
m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);

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
