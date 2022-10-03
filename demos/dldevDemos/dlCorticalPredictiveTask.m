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

TotalSize = [25, 50, 100, 200];
Currentsize = TotalSize(1);

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
ModelParametersPFC.NoiseRate = 10.5; % 15%
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
ModelParametersV4.NoiseRate = 14; % 20%
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
connection1.parameters={'gAMPA', .24, 'tauAMPA', 4, 'netcon', connectionWeigth1};

% supEV4->midPVPFC
connectionWeigth2 = 0.21*rand(dsCellV4.populations(1).size, dsCellPFC.populations(5).size) + 0.27;
connection2.direction = [dsCellV4.populations(1).name, '->', dsCellPFC.populations(5).name];

connection2.source = dsCellV4.populations(1).name;
connection2.target = dsCellPFC.populations(5).name;
connection2.mechanism_list={'iAMPActx'};
connection2.parameters={'gAMPA', .24, 'tauAMPA', 4, 'netcon', connectionWeigth2};

% deepEPFC->supSOMV4
connectionWeigth3 = 0.27*rand(dsCellPFC.populations(6).size, dsCellV4.populations(2).size) + 0.21;
connection3.direction = [dsCellPFC.populations(6).name, '->', dsCellV4.populations(2).name];

connection3.source = dsCellPFC.populations(6).name;
connection3.target = dsCellV4.populations(2).name;
connection3.mechanism_list={'iAMPActx'};
connection3.parameters={'gAMPA', .27, 'tauAMPA', 4, 'netcon', connectionWeigth3};

% deepEPFC->supEV4
set1 = [2:5, 16:19];
set2 = [6:9, 20:23];connectionWeigth4 = 0.27*rand(dsCellPFC.populations(6).size, dsCellV4.populations(1).size) + 0.21;
connection3.direction = [dsCellPFC.populations(6).name, '->', dsCellV4.populations(1).name];

connection3.source = dsCellPFC.populations(6).name;
connection3.target = dsCellV4.populations(1).name;
connection3.mechanism_list={'iAMPActx'};
connection3.parameters={'gAMPA', .27, 'tauAMPA', 4, 'netcon', connectionWeigth4};

%%% Finalization
dsModel = dlConnectModels({dsCellV4, dsCellPFC}, {connection1, connection2, connection3});

%% Create DynaLearn Class 
% Try to use this section only first time or If you have lost your file and
% you want a new model.

m = DynaLearn(dsModel, 'models/dlPredictiveCorticalCircuitModelLWK1', 'mex'); % ~10 min or less, MEXGEN or < 20 sec, RAWGEN.
m.dlSave(); % < 1sec

%% Load DynaLearn Class

m = DynaLearn(); % ~ 1sec
m = m.dlLoad('models/dlPredictiveCorticalCircuitModelLWK1'); % ~ 10sec, New larger model; keeping track of its activity in Gamma/Beta **

%% Trial: training  script preparation, 50-block and 50-trial

[trialParams1, trialParams2, trialParams3] = dlDemoThreePattern();

outputParams = [{'deepExPFC_V', 1:floor(ModelParametersPFC.NeDeep/3), [400 600] ...
    , 'alfp'}; {'deepExPFC_V',ceil(ModelParametersPFC.NeDeep/3):floor(2*ModelParametersPFC.NeDeep/3), ...
    [400 600], 'alfp'}; {'deepExPFC_V', ceil(2*ModelParametersPFC.NeDeep/3):ModelParametersPFC.NeDeep, [400 600], 'alfp'}; ...
    {'supExPFC_V', 1:ModelParametersPFC.NeSuperficial, [200 800], 'afr'}; ...
    {'midExPFC_V', 1:ModelParametersPFC.NeMid, [200 800], 'afr'}; ...
    {'deepExPFC_V', 1:ModelParametersPFC.NeDeep, [200 800], 'afr'}; ...
    {'supExV4_V', 1:ModelParametersV4.NeSuperficial, [200 800], 'afr'}; ...
    {'midExV4_V', 1:ModelParametersV4.NeMid, [200 800], 'afr'}; ...
    {'deepExV4_V', 1:ModelParametersV4.NeDeep, [200 800], 'afr'}];

targetParams1 = [{'TotalSpikesPenalty', 4:9, 160, 0.14}; {'MSE', 1, 40, 0.17}; {'MSE', 2, 20, 0.08}; {'MSE', 3, 20, 0.08}; {'Compare', [1, 2], 0, 0.04}; {'Compare', [1, 3], 0, 0.04}; {'Diff', [2, 3], 0, 0.01}]; % A 
targetParams2 = [{'TotalSpikesPenalty', 4:9, 160, 0.14}; {'MSE', 2, 40, 0.17}; {'MSE', 1, 20, 0.08}; {'MSE', 3, 20, 0.08}; {'Compare', [2, 1], 0, 0.04}; {'Compare', [2, 3], 0, 0.04}; {'Diff', [1, 3], 0, 0.01}]; % B
targetParams3 = [{'TotalSpikesPenalty', 4:9, 160, 0.14}; {'MSE', 3, 40, 0.17}; {'MSE', 2, 20, 0.08}; {'MSE', 1, 20, 0.08}; {'Compare', [3, 1], 0, 0.04}; {'Compare', [3, 2], 0, 0.04}; {'Diff', [1, 2], 0, 0.01}]; % C

dlInputParameters = {trialParams1, trialParams2, trialParams3};
dlTargetParameters = {targetParams1, targetParams2, targetParams3};
dlOutputParameters = outputParams;

TBdata = dlTrialBlockGenerator(dlInputParameters, dlTargetParameters, 50, 50);

dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
dlTrainOptions('dlEpochs') = 100; % % Number of epochs (A.K.A total iterations)
dlTrainOptions('dlBatchs') = 3; % If a scenario requires the training to be based on a group parameter (e.g mean of errors) use a dlBatch > 1 and set update mode later to batch. 
dlTrainOptions('dlLambda') = 1e-5; % Higher lambda means more changes based on error, lower may cause model to learn slower or nothing.
    
dlTrainOptions('dlCheckpoint') = 'true'; % If current step's error is higher based on a threshold, reload last optimal state and continue from that point
dlTrainOptions('dlCheckpointCoefficient') = 1.74; % A.K.A exploration rate
dlTrainOptions('dlUpdateMode') = 'batch'; % Update on each trial's result or based on batch group results

dlTrainOptions('dlLearningRule') = 'BioDeltaRule'; % Delta rule with a basic change based on biophysical properties 
dlTrainOptions('dlSimulationFlag') = 1; % If 0, will not run simulations (only for debugging purposes)
dlTrainOptions('dlOutputLogFlag') = 1; % If 0, will not keep outputs
dlTrainOptions('dlOfflineOutputGenerator') = 0; % If 1, will generate fake-random outputs (only for debugging purposes)

dlTrainOptions('dlAdaptiveLambda') = 0; % Adaptive lambda parameter; recommended for long simulations.
dlTrainOptions('dlLambdaCap') = 8e-3; % Only if Adaptive lambda is active, recommended to set a upper-bound (UB) or ignore to use default UB (0.01).

% Initial training on the model to reach a plausible local minimia like
% the task in the paper the model should also learn the basics of the task.
% We shortly train the model by cues to put it close to a local minimia.

dlTrainOptions('dlLambda') = 7e-4;
dlTrainOptions('dlEpochs') = 10;
dlTrainOptions('dlBatchs') = 3;

argsPowSpectRatio = struct();

argsPowSpectRatio.lf1 = 7;
argsPowSpectRatio.hf1 = 28;
argsPowSpectRatio.lf2 = 35;
argsPowSpectRatio.hf2 = 140;

dlTrainOptions('dlCustomLog') = "dlPowerSpectrumRatio"; % Name of a function which is in the path
dlTrainOptions('dlCustomLogArgs') = argsPowSpectRatio; % Arguments of your custom function

%% Pre-training : Train model to a "non-far" point.

clc;

dlTrainOptions('dlLambda') = 5e-5;
dlTrainOptions('dlEpochs') = 5;
m.dlOptimalError = 1e7;
% m.dlResetTraining();

tic;
m.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions); % <16 sec per trial
toc;

% disp(m.dlCurrentSessionValidTrials);
% disp(size(m.dlErrorsLog));

%% Block-trial phase

% TODO: Remove "far" trials from training by backtrack.
% TODO: Add raster plotter.
% TODO: Add live runs (continouos task) option.

clc;
dlTrainOptions('dlLambda') = 1e-7;
dlTrainOptions('dlAdaptiveLambda') = 0; % Adaptive lambda parameter; recommended for long simulations.
dlTrainOptions('dlUpdateMode') = 'trial';
dlTrainOptions('dlEpochs') = 10;
dlTrainOptions('dlBatchs') = 50;

m.dlResetTraining();
argsPSR = struct();

argsPSR.lf1 = 7;
argsPSR.hf1 = 28;
argsPSR.lf2 = 35;
argsPSR.hf2 = 140;

dlTrainOptions('dlCustomLog') = "dlPowerSpectrumRatio"; % Name of a function which is in the path
dlTrainOptions('dlCustomLogArgs') = argsPSR; % Arguments of your custom function

disp("----------A-----------");
dlTrainOptions('dlCheckpointCoefficient') = 1.7;
m.dlTrain(TBdata.B1, dlOutputParameters, TBdata.T1, dlTrainOptions);

disp("----------U1----------");
m.dlErrorsLog = [m.dlErrorsLog, -1];  
m.dlOptimalError = 1e9;
dlTrainOptions('dlCheckpointCoefficient') = 4.7;
m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);

disp("----------B-----------");
m.dlErrorsLog = [m.dlErrorsLog, -1]; 
m.dlOptimalError = 1e9;
dlTrainOptions('dlCheckpointCoefficient') = 1.7; 
m.dlTrain(TBdata.B2, dlOutputParameters, TBdata.T2, dlTrainOptions);

disp("----------U2----------");
m.dlErrorsLog = [m.dlErrorsLog, -1]; 
m.dlOptimalError = 1e9;
dlTrainOptions('dlCheckpointCoefficient') = 4.7; 
m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);

disp("----------C----------");
m.dlErrorsLog = [m.dlErrorsLog, -1]; 
m.dlOptimalError = 1e9;
dlTrainOptions('dlCheckpointCoefficient') = 1.7;
m.dlTrain(TBdata.B3, dlOutputParameters, TBdata.T3, dlTrainOptions);

disp("----------U3----------");
m.dlErrorsLog = [m.dlErrorsLog, -1]; 
m.dlOptimalError = 1e9;
dlTrainOptions('dlCheckpointCoefficient') = 4.7;
m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);

%% Errors log plot

% clc;

wk = 2;
figure('position', [0, 0, 1400, 700]);
n = max(size(m.dlErrorsLog));
x = zeros(1, ceil(n/wk));

for i = 0:wk:n-wk
    x(ceil((i+1)/wk)) = min(m.dlErrorsLog(i+1:i+wk));
end
tlab = ["A", "U1", "B", "U2", "C", "U3"];
w = 100;
errorcap = max(x);

for k = 1:6
    
    for l = 0:5
        fill([l*w, l*w, l*w+w, l*w+w], [0, errorcap*1.2, errorcap*1.2, 0], [sin(l*0.2), 1, cos(l*0.2)]);hold('on');
        text(l*w+10, errorcap*1.1, tlab(l+1));
        ylim([0 errorcap*1.2]);
    end

end

xlim([0 600]);
plot(x);grid("on");
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

w = 100;
cnt = 1;
tlab = ["A", "U1", "B", "U2", "C", "U3"];
f = figure("Position", [0 0 1500 1000]);

set1 = [2:5, 16:19];
set2 = [6:9, 20:23];

for i = 1:4

    for l = 0:5    

        subplot(2, 2, i);
        fill([l*w, l*w, l*w+w, l*w+w], [5, 40, 40, 5], [sin(l*0.2), 1, cos(l*0.2)], 'HandleVisibility','off');
        hold('on');

    end

end

for i = set2

    q = m.dlCustomLog(i, 601:1200);  
    subplot(2, 2, mod(cnt-1, 4)+1);
    plot(q, 'DisplayName', m.dlCustomLogLabel{i});
    legend("Location", "southwest"); 

    if i > 15

        ql = min(min(q), min(m.dlCustomLog(i-14, :)));
        qu = max(max(q), max(m.dlCustomLog(i-14, :)));
        ylim([ql-1 qu+1]);

        for l = 0:5
            text((l+0.5)*w, qu+0.5, tlab(l+1));
        end

        xlabel(m.dlCustomLogLabel{i} + " and " + m.dlCustomLogLabel{i-14});
    
    end

    hold("on");
    grid("on");
    cnt = cnt + 1;

end

sgtitle("Relative Beta to Gamma power band ratio.");

%%

clc;

w = 100;
cnt = 1;
tlab = ["A", "U1", "B", "U2", "C", "U3"];
figure("Position", [0 0 1500 1000]);
% f.Name = "aaaaaaa";

for i = 1:4

    for l = 0:5    

        subplot(2, 2, i);
        fill([l*w, l*w, l*w+w, l*w+w], [5, 30, 30, 5], [sin(l*0.2), 1, cos(l*0.2)], 'HandleVisibility','off');
        hold('on');

    end

end

for i = [6:9, 20:23]

    q = m.dlCustomLog(i, :);  
    subplot(2, 2, mod(cnt-1, 4)+1);
    plot(q, 'DisplayName', m.dlCustomLogLabel{i});
    legend("Location", "southwest"); 

    if i > 15

        ql = min(min(q, m.dlCustomLog(i-14, :)));
        qu = max(max(q, m.dlCustomLog(i-14, :)));
        ylim([ql-1 qu+1]);
                
        for l = 0:5
            text((l+0.5)*w, qu+0.5, tlab(l+1));
        end

        xlabel(m.dlCustomLogLabel{i} + " and " + m.dlCustomLogLabel{i-14});
    
    end

    hold("on");
    grid("on");
    cnt = cnt + 1;

end

sgtitle("Relative Beta to Gamma power band ratio.");

%% End