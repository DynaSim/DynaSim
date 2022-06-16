%% Simulation of predictive task implemented on DynaSim/DynaLearn via Reinforcement learning 

% Based on [Bastos2020:Layer and rhythm specificity for predictive routing, A.M Bastos et. al, 2020]
% Requirement: DynaSim/dev latest version
% Minimum MatLab version: R2019a
% 5th of July, 2022
% @HNXJ

%% Model parameters

clear;clc;
Ne = 24;Ni = 6;Nin = 6;NoiseRate = 5;
s3 = dlModelPredictivePFC(Ne, Ni, Nin, NoiseRate); % Predictive PFC model with specific parameters

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

TBdata = dlTrialBlockGenerator(dlInputParameters, dlTargetParameters, 50, 50);

dlTrainOptions = containers.Map();
dlTrainOptions('dlEpochs') = 1;
dlTrainOptions('dlBatchs') = 50;
dlTrainOptions('dlLambda') = 1e-5;
    
dlTrainOptions('dlCheckpoint') = 'true';
dlTrainOptions('dlCheckpointCoefficient') = 2.74; % A.K.A exploration rate 
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

dlTrainOptions('dlEpochs') = 97;
dlTrainOptions('dlBatchs') = 3;
m.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

%% Block-trial phase

clc;
m.dlTrain(B1, dlOutputParameters, T1, dlTrainOptions);
m.dlTrain(TrB, dlOutputParameters, TrT, dlTrainOptions);
m.dlTrain(B2, dlOutputParameters, T2, dlTrainOptions);
m.dlTrain(TrB, dlOutputParameters, TrT, dlTrainOptions);
m.dlTrain(B3, dlOutputParameters, T3, dlTrainOptions);
m.dlTrain(TrB, dlOutputParameters, TrT, dlTrainOptions);

%% Errors log plot

% clc;
m.dlPlotBatchErrors(1);

%% Plot Local-field potentials

clc;
m.dlPlotAllPotentials('lfp');

%% Run a simulation (without training)

for i = 1:3
    m.dlRunSimulation(dlInputParameters{i}, dlOutputParameters);
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
