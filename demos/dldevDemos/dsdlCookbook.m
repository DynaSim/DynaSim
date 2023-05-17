%% DynaSim/DynaLearn cookbook
% github@hnxj; github@dynasim

% Run this after specifying your path to dynasim.
clear;close all;clc;
PathToDynaSim = 'D:\Works\Computational'; % Change it based on your local path for dynasim
cd(PathToDynaSim);
addpath(genpath('DynaSim'));
cd('DynaSim');

%% Single population with signle cell optimaztion examples
% Call function to make a dynasim model. Modify if you want more neurons or different noise level.

% DynaSim struct
populationName = "SingleHH";
ModelParameters = struct();
ModelParameters.Size = 1; % Cell count
ModelParameters.Noise = 7; % in mV, 7 is about 10% variance.
dsM = dlSinglePopulationNS(ModelParameters, populationName);

%% Model details
% Population names and connections. Consider that in this model, each cell has soma and dendrite.
dsM.populations.name
dsM.connections.direction

%% Passing to DynaLearn
% Assign an ID number if you want to avoid conflicts.

id = 0;
dlM = DynaLearn(dsM, char("dlModels/dlHHNS" + num2str(id)), 'mex', populationName);
dlM.dlSave();

%% Load (ONLY IF ALREADY EXISTS)

id = 0;
dlM = DynaLearn.dlLoader(char("dlModels/dlHHNS" + num2str(id)));

%% Define optimization parameters

inputParams = dlNullInputs(1000); % Simulation duration for each trial

outputParams = [{['ESOMA_SingleHH', '_V'], 1:1, [100 1000], 'av'}; ...
    {['EDEND_SingleHH', '_V'], 1:1, [100 1000], 'av'}; ...
    {['ESOMA_SingleHH', '_V'], 1:1, [100 1000], 'afr'}; ...
    {['EDEND_SingleHH', '_V'], 1:1, [100 1000], 'afr'}];

targetParams = [{'MSE', 1, -57, .4, 0, 0, 0, 0}; ...
    {'MSE', 2, -57, .4, 0, 0, 0, 0}; ...
    {'Compare', [1, 2], 0, .2, 0, 0, 0, 0}; ...
    {'MSE', 3, 10, .4, 0, 0, 0, 0}; ...
    {'MSE', 4, 10, .4, 0, 0, 0, 0}; ...
    {'Compare', [3, 4], 0, .2, 0, 0, 0, 0}];

dlInputParameters = {inputParams};
dlTargetParameters = {targetParams};
dlOutputParameters = outputParams;

dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
dlTrainOptions('dlEpochs') = 100; % % Number of epochs (A.K.A total iterations)
dlTrainOptions('dlBatchs') = 1; % If a scenario requires the training to be based on a group parameter (e.g mean of errors) use a dlBatch > 1 and set update mode later to batch. 
dlTrainOptions('dlLambda') = 4e-7; % Higher lambda means more changes based on error, lower may cause model to learn slower or nothing.
    
dlTrainOptions('dlCheckpoint') = 'true'; % If current step's error is higher based on a threshold, reload last optimal state and continue from that point
dlTrainOptions('dlCheckpointCoefficient') = 4.047; % A.K.A exploration rate
dlTrainOptions('dlCheckpointLengthCap') = 11; % If more than 7 steps with no progress passed, return to last checkpoint.
dlTrainOptions('dlUpdateMode') = 'trial'; % Update on each trial's result or based on batch group results

% dlTrainOptions('dlLearningRule') = 'BioDeltaRule'; % Delta rule with a basic change based on biophysical properties 
dlTrainOptions('dlSimulationFlag') = 1; % If 0, will not run simulations (only for debugging purposes)
dlTrainOptions('dlOutputLogFlag') = 0; % If 0, will not keep outputs
dlTrainOptions('dlOfflineOutputGenerator') = 0; % If 1, will generate fake-random outputs (only for debugging purposes)

dlTrainOptions('dlAdaptiveLambda') = 1; % Adaptive lambda parameter; recommended for long simulations.
dlTrainOptions('dlLambdaCap') = 1.1; % Only if Adaptive lambda is active, recommended to set a upper-bound (UB) or ignore to use default UB (0.01).
dlTrainOptions('dlExcludeDiverge') = 1; % Exclude non-optimals from model log
% dlTrainOptions('dlTrainExcludeList') = {'Stim', ['DeepISOM', ModelName, '->'], ['DeepIPV', ModelName, '->'], ['MidIPV', ModelName, '->'], ['supIPV', ModelName, '->'], ['supISOM', ModelName, '->']}; % Exclude populations from training

argsPowSpectRatio1 = struct();
argsPowSpectRatio2 = struct();
argsPowSpectRatio3 = struct();
argsPowSpectRatio4 = struct();

argsPowSpectRatio1.lf1 = 2;
argsPowSpectRatio1.hf1 = 6;
argsPowSpectRatio2.lf1 = 8;
argsPowSpectRatio2.hf1 = 14;

argsPowSpectRatio3.lf1 = 15;
argsPowSpectRatio3.hf1 = 30;
argsPowSpectRatio4.lf1 = 40;
argsPowSpectRatio4.hf1 = 90;

dlTrainOptions('dlCheckpointLengthCap') = 10;
dlTrainOptions('dlTrainExcludeList') = {'Stim'};
dlTrainOptions('dlCustomLog') = ["dlEPowerSpectrum", "dlEPowerSpectrum", "dlEPowerSpectrum", "dlEPowerSpectrum", "dlLFPaxLog"]; % Name of a function which is in the path
dlTrainOptions('dlCustomLogArgs') = [argsPowSpectRatio1, argsPowSpectRatio2, argsPowSpectRatio3, argsPowSpectRatio4, argsPowSpectRatio1, argsPowSpectRatio1]; % Arguments of your custom function

%% Optimization based on cell connection (dendrite-soma interaction)

dlTrainOptions('dlLambda') = 1e-9;
dlTrainOptions('dlAdaptiveLambda') = 1; % Adaptive lambda parameter; recommended for long simulations.
dlTrainOptions('dlUpdateMode') = 'trial';
dlTrainOptions('dlLearningRule') = 'EnhancedDeltaRule';

dlTrainOptions('dlEpochs') = 10;
dlTrainOptions('dlBatchs') = 1;
dlTrainOptions('dlEnhancedMomentum') = 0.2;
dlTrainOptions('dlCheckpointCoefficient') = 2.074;

dlM.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 1;
dlM.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);
dlM.dlSave();

%% Optimization based on all variables (dendrite-soma interaction plus conductance(g), time constant(tau) and ...)

dlTrainOptions('dlLambda') = 1e-9;
dlTrainOptions('dlAdaptiveLambda') = 1; % Adaptive lambda parameter; recommended for long simulations.
dlTrainOptions('dlUpdateMode') = 'trial';
dlTrainOptions('dlLearningRule') = 'GeneralEnhancedDeltaRule';

dlTrainOptions('dlEpochs') = 10000;
dlTrainOptions('dlBatchs') = 1;
dlTrainOptions('dlEnhancedMomentum') = 0.2;
dlTrainOptions('dlCheckpointCoefficient') = 2.074;

dlM.dlLoadOptimal();
dlM.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 1;
dlM.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

%% END  

