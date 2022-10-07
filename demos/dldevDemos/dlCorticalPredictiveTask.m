%% Dual-Area simulation of predictive task implemented on DynaSim/DynaLearn via Reinforcement learning 

% Based on papers: 
% [Layer and rhythm specificity for predictive routing, A.M Bastos et. al, 2020]
% [Top-Down beta rhythms support selective attention via interlaminar interaction: a model, J.H Lee et. al, 2013]
% [DynaSim: a MATLAB toolbox for neural modeling and simulation, J.S Sherfey et. al, 2018]

% Requirement: DynaSim/dev_2022 latest version
% Minimum MatLab version: R2021a
% @Septembre2022
% @HNXJ

%% Model parameters

clear;clc;

model_size_id = 1;
TotalSize = [50, 50, 50, 50];
Currentsize = TotalSize(model_size_id);

%%% Create model parameters struct
ModelParametersPFC = struct();

%%% Area PFC layer sizes (relative)
ModelParametersPFC.NeSuperficial = ceil(0.3*Currentsize);
ModelParametersPFC.NSomSuperficial = ceil(0.02*Currentsize);
ModelParametersPFC.NPvSuperficial = ceil(0.04*Currentsize);
ModelParametersPFC.NeMid = ceil(0.14*Currentsize);
ModelParametersPFC.NSomMid = 0;
ModelParametersPFC.NPvMid = ceil(0.02*Currentsize);
ModelParametersPFC.NeDeep = ceil(0.45*Currentsize);
ModelParametersPFC.NSomDeep = ceil(0.02*Currentsize);
ModelParametersPFC.NPvDeep = ceil(0.01*Currentsize);

ModelParametersPFC.Nin = 6;
ModelParametersPFC.Nout = 6;
ModelParametersPFC.NoiseRate = 4; % 6%
ModelParametersPFC.Nstim = 3;

%%% Area V4 layer sizes (relative) 
ModelParametersV4.NeSuperficial = ceil(0.25*Currentsize);
ModelParametersV4.NSomSuperficial = ceil(0.03*Currentsize);
ModelParametersV4.NPvSuperficial = ceil(0.07*Currentsize);
ModelParametersV4.NeMid = ceil(0.12*Currentsize);
ModelParametersV4.NSomMid = 0;
ModelParametersV4.NPvMid = ceil(0.03*Currentsize);
ModelParametersV4.NeDeep = ceil(0.45*Currentsize);
ModelParametersV4.NSomDeep = ceil(0.03*Currentsize);
ModelParametersV4.NPvDeep = ceil(0.02*Currentsize);

ModelParametersV4.Nin = 6;
ModelParametersV4.Nout = 6;
ModelParametersV4.NoiseRate = 7; % 10%
ModelParametersV4.Nstim = 3;

%%% Call Laminar Cortex Constructor Functions
dsCellPFC = dlLaminarCortexNetLWK(ModelParametersPFC, 'PFC'); % Laminar PFC model with specific parameters
dsCellV4 = dlLaminarCortexNetLWK(ModelParametersV4, 'V4'); % Laminar V4 model with specific parameters

%%% Connecting two models: Define connections between two areas
% supEV4->midEPFC
connectionWeigth1 = 0.21*rand(dsCellV4.populations(1).size, dsCellPFC.populations(4).size) + 0.27;
connection1.direction = [dsCellV4.populations(1).name, '->', dsCellPFC.populations(4).name];

connection1.source = dsCellV4.populations(1).name;
connection1.target = dsCellPFC.populations(4).name;
connection1.mechanism_list={'iAMPActx'};
connection1.parameters={'gAMPA', .3, 'tauAMPA', 1, 'netcon', connectionWeigth1};

% supEV4->midPVPFC
connectionWeigth2 = 0.21*rand(dsCellV4.populations(1).size, dsCellPFC.populations(5).size) + 0.27;
connection2.direction = [dsCellV4.populations(1).name, '->', dsCellPFC.populations(5).name];

connection2.source = dsCellV4.populations(1).name;
connection2.target = dsCellPFC.populations(5).name;
connection2.mechanism_list={'iAMPActx'};
connection2.parameters={'gAMPA', .3, 'tauAMPA', 1, 'netcon', connectionWeigth2};

% deepEPFC->supSOMV4
connectionWeigth3 = 0.27*rand(dsCellPFC.populations(6).size, dsCellV4.populations(2).size) + 0.21;
connection3.direction = [dsCellPFC.populations(6).name, '->', dsCellV4.populations(2).name];

connection3.source = dsCellPFC.populations(6).name;
connection3.target = dsCellV4.populations(2).name;
connection3.mechanism_list={'iAMPActx'};
connection3.parameters={'gAMPA', .9, 'tauAMPA', 1, 'netcon', connectionWeigth3};

% deepEPFC->supEV4
connectionWeigth4 = 0.27*rand(dsCellPFC.populations(6).size, dsCellV4.populations(1).size) + 0.21;
connection4.direction = [dsCellPFC.populations(6).name, '->', dsCellV4.populations(1).name];

connection4.source = dsCellPFC.populations(6).name;
connection4.target = dsCellV4.populations(1).name;
connection4.mechanism_list={'iAMPActx'};
connection4.parameters={'gAMPA', .9, 'tauAMPA', 1, 'netcon', connectionWeigth4};

%%% Finalization
dsModel = dlConnectModels({dsCellV4, dsCellPFC}, {connection1, connection2, connection3, connection4});

%% Create DynaLearn Class 
% Try to use this section only first time or If you have lost your file and
% you want a new model.

m = DynaLearn(dsModel, char("models/dlPredictiveCorticalCircuitModelLWK" + string(model_size_id)), 'mex'); % ~10 min or less, MEXGEN or < 20 sec, RAWGEN.
m.dlSave(); % < 1sec

%% Load DynaLearn Class

clc;
% model_size_id = 1;
m = DynaLearn(); % ~ 1sec
m = m.dlLoad(char("models/dlPredictiveCorticalCircuitModelLWK" + string(model_size_id))); % ~ 10sec, New larger model; keeping track of its activity in Gamma/Beta **

%% Trial: training  script preparation, 50-block and 50-trial

clc;

[trialParams1, trialParams2, trialParams3] = dlDemoThreePattern();

outputParams = [{'deepExPFC_V', 1:floor(ModelParametersPFC.NeDeep/3), [400 750] ...
    , 'afr'}; {'deepExPFC_V',ceil(ModelParametersPFC.NeDeep/3):floor(2*ModelParametersPFC.NeDeep/3), ...
    [400 750], 'afr'}; {'deepExPFC_V', ceil(2*ModelParametersPFC.NeDeep/3):ModelParametersPFC.NeDeep, [400 750], 'afr'}; ...
    {'supExPFC_V', 1:ModelParametersPFC.NeSuperficial, [300 800], 'afr'}; ...
    {'midExPFC_V', 1:ModelParametersPFC.NeMid, [300 800], 'afr'}; ...
    {'deepExPFC_V', 1:ModelParametersPFC.NeDeep, [300 800], 'afr'}; ...
    {'supExV4_V', 1:ModelParametersV4.NeSuperficial, [300 800], 'afr'}; ...
    {'midExV4_V', 1:ModelParametersV4.NeMid, [300 800], 'afr'}; ...
    {'deepExV4_V', 1:ModelParametersV4.NeDeep, [300 800], 'afr'}];

targetParams1 = [{'TotalSpikesPenalty', 4:7, 100, 0.4}; {'MSE', 1, 50, 0.05}; {'MSE', 2, 25, 0.01}; {'MSE', 3, 25, 0.01}; {'Compare', [1, 2], 0, 0.2}; {'Compare', [1, 3], 0, 0.2}; {'Diff', [2, 3], 0, 0.01}]; % A 
targetParams2 = [{'TotalSpikesPenalty', 4:7, 100, 0.4}; {'MSE', 2, 50, 0.05}; {'MSE', 1, 25, 0.01}; {'MSE', 3, 25, 0.01}; {'Compare', [2, 1], 0, 0.2}; {'Compare', [2, 3], 0, 0.2}; {'Diff', [1, 3], 0, 0.01}]; % B
targetParams3 = [{'TotalSpikesPenalty', 4:7, 100, 0.4}; {'MSE', 3, 50, 0.05}; {'MSE', 2, 25, 0.01}; {'MSE', 1, 25, 0.01}; {'Compare', [3, 2], 0, 0.2}; {'Compare', [3, 1], 0, 0.2}; {'Diff', [2, 1], 0, 0.01}]; % C

dlInputParameters = {trialParams1, trialParams2, trialParams3};
dlTargetParameters = {targetParams1, targetParams2, targetParams3};
dlOutputParameters = outputParams;

TBdata = dlTrialBlockGenerator(dlInputParameters, dlTargetParameters, 25, 25);

dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
dlTrainOptions('dlEpochs') = 100; % % Number of epochs (A.K.A total iterations)
dlTrainOptions('dlBatchs') = 3; % If a scenario requires the training to be based on a group parameter (e.g mean of errors) use a dlBatch > 1 and set update mode later to batch. 
dlTrainOptions('dlLambda') = 1e-5; % Higher lambda means more changes based on error, lower may cause model to learn slower or nothing.
    
dlTrainOptions('dlCheckpoint') = 'true'; % If current step's error is higher based on a threshold, reload last optimal state and continue from that point
dlTrainOptions('dlCheckpointCoefficient') = 2.047; % A.K.A exploration rate
dlTrainOptions('dlCheckpointLengthCap') = 7; % If more than 7 steps with no progress passed, return to last checkpoint.
dlTrainOptions('dlUpdateMode') = 'batch'; % Update on each trial's result or based on batch group results

dlTrainOptions('dlLearningRule') = 'BioDeltaRule'; % Delta rule with a basic change based on biophysical properties 
dlTrainOptions('dlSimulationFlag') = 1; % If 0, will not run simulations (only for debugging purposes)
dlTrainOptions('dlOutputLogFlag') = 1; % If 0, will not keep outputs
dlTrainOptions('dlOfflineOutputGenerator') = 0; % If 1, will generate fake-random outputs (only for debugging purposes)

dlTrainOptions('dlAdaptiveLambda') = 1; % Adaptive lambda parameter; recommended for long simulations.
dlTrainOptions('dlLambdaCap') = 1.1; % Only if Adaptive lambda is active, recommended to set a upper-bound (UB) or ignore to use default UB (0.01).
dlTrainOptions('dlExcludeDiverge') = 1; % Exclude non-optimals from model log
dlTrainOptions('dlTrainExcludeList') = {'Stim'}; % Exclude populations from training

dlTrainOptions('dlLambda') = 7e-4;
dlTrainOptions('dlEpochs') = 10;
dlTrainOptions('dlBatchs') = 3;

argsPowSpectRatio = struct();
argsNull = [];

argsPowSpectRatio.lf1 = 7;
argsPowSpectRatio.hf1 = 28;
argsPowSpectRatio.lf2 = 42;
argsPowSpectRatio.hf2 = 84;

dlTrainOptions('dlCustomLog') = ["dlPowerSpectrumRatio", "dlAccuracyBastos2020Task"]; % Name of a function which is in the path
dlTrainOptions('dlCustomLogArgs') = [argsPowSpectRatio, argsNull]; % Arguments of your custom function

%% Pre-training : Train model to a "non-far" point.

clc;

dlTrainOptions('dlLambda') = 4e-2;
dlTrainOptions('dlEpochs') = 5;
dlTrainOptions('dlCheckpointCoefficient') = 1.4; 
dlTrainOptions('dlCheckpointLengthCap') = 14;

m.dlOptimalError = 1e7;
m.dlResetTraining();

tic;
% m.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions); % <16 sec per trial
toc;


%% Block-trial phase

% TODO: Add average U2B/B2U
% TODO: 

% TODO: Add raster plotter.
% TODO: Add live runs (continouos task) option.

clc;

dlTrainOptions('dlLambda') = 1e-9; % 1e-11(1) -> 1e-4 (4)
dlTrainOptions('dlAdaptiveLambda') = 0; % Adaptive lambda parameter; recommended for long simulations.
dlTrainOptions('dlUpdateMode') = 'trial';
dlTrainOptions('dlLearningRule') = 'BioDeltaRule';

dlTrainOptions('dlTrainExcludeList') = {'Stimuli'};
dlTrainOptions('dlCheckpointLengthCap') = 15;
dlTrainOptions('dlEpochs') = 2;
dlTrainOptions('dlBatchs') = 25;

dlTrainOptions('dlExcludeDiverge') = 1;
CheckCoeff = 1.25;

m.dlResetTraining();
argsPSR = struct();

argsPSR.lf1 = 7;
argsPSR.hf1 = 28;
argsPSR.lf2 = 35;
argsPSR.hf2 = 140;

dlTrainOptions('dlCustomLog') = ["dlPowerSpectrumRatio", "dlAccuracyBastos2020Task"]; % Name of a function which is in the path
dlTrainOptions('dlCustomLogArgs') = [argsPSR, argsNull]; % Arguments of your custom function

for cnt = 1:1

    disp("----------A-----------");
    dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff;
    m.dlTrain(TBdata.B1, dlOutputParameters, TBdata.T1, dlTrainOptions);
    
    disp("----------U1----------");
    m.dlErrorsLog = [m.dlErrorsLog, -1];  
    m.dlOptimalError = 1e9;
    dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff*1.2;
    m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);
    
    disp("----------B-----------");
    m.dlErrorsLog = [m.dlErrorsLog, -1]; 
    m.dlOptimalError = 1e9;
    dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff; 
    m.dlTrain(TBdata.B2, dlOutputParameters, TBdata.T2, dlTrainOptions);
    
    disp("----------U2----------");
    m.dlErrorsLog = [m.dlErrorsLog, -1]; 
    m.dlOptimalError = 1e9;
    dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff*1.2; 
    m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);
    
    disp("----------C----------");
    m.dlErrorsLog = [m.dlErrorsLog, -1]; 
    m.dlOptimalError = 1e9;
    dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff;
    m.dlTrain(TBdata.B3, dlOutputParameters, TBdata.T3, dlTrainOptions);
    
    disp("----------U3----------");
    m.dlErrorsLog = [m.dlErrorsLog, -1]; 
    m.dlOptimalError = 1e9;
    dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff*1.2;
    m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);

end

%% Errors log plot

% clc;

wk = 1;
figure('position', [0, 0, 1400, 700]);
n = max(size(m.dlErrorsLog));
x = zeros(1, ceil(n/wk));

for i = 0:wk:n-wk
    x(ceil((i+1)/wk)) = min(m.dlErrorsLog(i+1:i+wk));
end

tlab = ["A", "U1", "B", "U2", "C", "U3", "A", "U4", "B", "U5", "C", "U6"];
w = 50;
errorcap = max(x);

for k = 1:6
    
    for l = 0:11

        fill([l*w, l*w, l*w+w, l*w+w], [0, errorcap*1.2, errorcap*1.2, 0], [sin(l*0.1), 1, cos(l*0.1)]);hold('on');
        text(l*w+10, errorcap*1.1, tlab(l+1));
        ylim([0 errorcap*1.2]);

    end

end

xlim([0 600]);
plot(x(x>0));grid("on");
title("Errors in trials");

%% Plot Local-field potentials

clc;
m.dlPlotAllPotentials('lfp');

%% Raw run of a simulation (without training)

clc;
tic;
for i = 1:3

    m.dlRunSimulation(dlInputParameters{i}, dlOutputParameters);

end
toc;

%% Beta to Gamma ratio (relative) in trials during predictive task training

clc;

w = 200;
cnt = 1;
tlab = ["A", "U1", "B", "U2", "C", "U3", "A", "U4", "B", "U5", "C", "U6"];
f = figure("Position", [0 0 1500 1000]);

set1 = [2:5, 16:19];
set2 = [6:9, 20:23];

ratioLog = {m.dlCustomLog(1, :)};
ratioLog = cell2mat(ratioLog{1});

for i = 1:4

    for l = 0:11

        subplot(2, 2, i);
        fill([l*w, l*w, l*w+w, l*w+w], [5, 40, 40, 5], [sin(l*0.1), 1, cos(l*0.1)], 'HandleVisibility','off');
        hold('on');

    end

end

dispLabels = m.dlCustomLogLabel{1};

for i = set2

    subplot(2, 2, mod(cnt-1, 4)+1);
    plot(ratioLog(i, :), 'DisplayName', dispLabels{i});
    legend("Location", "southwest");

    if i > 15

        q = ratioLog(i-14, :);
        q2 = ratioLog(i, :);

        ql = min(min(q), min(q2));
        qu = max(max(q), max(q2));
        ylim([ql-1 qu+1]);

        for l = 0:11
            text((l+0.5)*w, qu+0.5, tlab(l+1));
        end

        xlabel(dispLabels{i} + " and " + dispLabels{i-14});
    
    end

    hold("on");
    grid("on");
    cnt = cnt + 1;

end

sgtitle("Relative Beta to Gamma power band ratio.");

%% AVG Beta to Gamma ratio (relative) in trials during predictive task training

clc;

w = 100;
cnt = 1;
tlab = ["Predictable", "Unpredictable"];
f = figure("Position", [0 0 1500 1000]);

set1 = [2:5, 16:19];
set2 = [6:9, 20:23];

ratioLog = {m.dlCustomLog(1, :)};
ratioLog = cell2mat(ratioLog{1});

for i = 1:4

    for l = 0:1

        subplot(2, 2, i);
        fill([l*w, l*w, l*w+w, l*w+w], [0, 40, 40, 0], [sin(l*0.9), 0.75, cos(l*0.001)], 'HandleVisibility','off');
        hold('on');

    end

end

dispLabels = m.dlCustomLogLabel{1};
ratioLogAvg = zeros(size(ratioLog, 1), 200);
K = 11;

for j = 1:K

    ratioLogAvg(:, :) = ratioLogAvg + ratioLog(:, (j-1)*200 + 1:j*200)/K;

end

for j = 1:size(ratioLogAvg, 1)

    ratioLogAvg(j, 1:196) = conv(ratioLogAvg(j, :), ones(1, 5)/5, "valid");

end

for i = set2



    subplot(2, 2, mod(cnt-1, 4)+1);
    plot(ratioLogAvg(i, :), 'DisplayName', dispLabels{i});
    legend("Location", "southwest");

    if i > 15

        q = ratioLogAvg(i-14, :);
        q2 = ratioLogAvg(i, :);

        ql = min(min(q), min(q2));
        qu = max(max(q), max(q2));
        ylim([ql-1 qu+1]);

        for l = 0:1
            text((l+0.5)*w, qu+0.5, tlab(l+1));
        end

        title(dispLabels{i} + " and " + dispLabels{i-14});
    
    end

    hold("on");
    grid("on");
    cnt = cnt + 1;

end

sgtitle("Average Beta/Gamma power band ratio across (14) different trial transitions");

%% End