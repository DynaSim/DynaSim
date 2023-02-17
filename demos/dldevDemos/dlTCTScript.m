%% Single-Area simulation of predictive task implemented on DynaSim/DynaLearn via Reinforcement learning 

% Based on papers: 
% [Layer and rhythm specificity for predictive routing, A.M Bastos et. al, 2020]
% [Top-Down beta rhythms support selective attention via interlaminar interaction: a model, J.H Lee et. al, 2013]
% [DynaSim: a MATLAB toolbox for neural modeling and simulation, J.S Sherfey et. al, 2018]

% Requirement: DynaSim/dev_2022 latest version
% Minimum MatLab version: R2021a
% @Septembre2022
% @HNXJ

%% Debugging box

clear;clc;
cd('/home/hamednejat/Works/dsm');
addpath(genpath('DynaSim'));
cd('DynaSim');

%% Debugging area

ModelParametersPFC = struct();
PFCSize = 64;

%%% Area PFC layer sizes (relative)
ModelParametersPFC.NeSuperficial = ceil(0.24*PFCSize);
ModelParametersPFC.NPvSuperficial = ceil(0.02*PFCSize);
ModelParametersPFC.NSomSuperficial = ceil(0.02*PFCSize);

ModelParametersPFC.NeMid = ceil(0.12*PFCSize);
ModelParametersPFC.NSomMid = ceil(0.00*PFCSize);
ModelParametersPFC.NPvMid = ceil(0.04*PFCSize);

ModelParametersPFC.NeDeep = ceil(0.34*PFCSize);
ModelParametersPFC.NSomDeep = ceil(0.05*PFCSize);
ModelParametersPFC.NPvDeep = ceil(0.05*PFCSize);

ModelParametersPFC.Nin = 6;
ModelParametersPFC.Nout = 6;
ModelParametersPFC.NoiseRate = 5.7;
ModelParametersPFC.Nstim = 3;

%%% Call Laminar Cortex Constructor Functions
dsCellPFC = dlLaminarCortexNetLWK(ModelParametersPFC, 'PFC'); % Laminar PFC model with specific parameters

%%% Finalization
dsModel = dsCellPFC;
simulator_options = {'tspan', [0 750], 'solver', 'rk1', 'dt', .01, ...
            'downsample_factor', 10, 'verbose_flag', 1};
dsData = dsSimulate(dsModel, simulator_options{:});
dsPlot(dsData);

%% Plot sample response

clear;clc;
id = 1;
m = DynaLearn(); % ~ 1sec
m = m.dlLoad(char("dlModels/dlTCTModels" + string(id))); % ~ 10sec, New larger model; keeping track of its activity in Gamma/Beta **

trialParamsTemp = containers.Map();
trialParamsTemp('tspan') = [0 1500];

m.dlUpdateParams(trialParamsTemp);
m.dlSimulate();
m.dlPlotAllPotentials('raster');
% opts = struct();opts.name = "test1";
% m.dlPlotAllPotentials('lfpsave', opts);   

%%

figure('Position', [0, 0, 1700, 1400]);

for i = 5

    load("dlModels/dlTCTModels2/solve/signals_test" + num2str(i) + ".mat");

%     subplot(3, 3, i);
    imagesc(imgRaster>4);
    colormap("gray");

end

%% AutoRunSc (SinglePFC)

clear;clc;
 
TotalSize = ones(1, 10)*128;
noise_rate = 4.7;

for model_size_id = 1:1

    Currentsize = TotalSize(model_size_id);
    m = dlTCTPerformerDynamics(1, Currentsize, model_size_id, noise_rate);

end

%% Continue learning

clc;

currentSize = 128;
reGenerate = 0;
noiseRate = 5;
modelID = 1;

for i = 1:10

    m = dlTCTPerformerDynamics(reGenerate, currentSize, modelID, noiseRate);
    opts = struct();opts.name = "test" + num2str(i);
    m.dlPlotAllPotentials('lfpsave', opts);

end

%%

x = m.dlErrorsLog;
x = smooth(x, 40);
x(end:end+1000) = m.dlErrorsLog(end-1000:end);
figure("Position", [0 0 1500 1000]);
plot(x);
title("Error/Loss (MSE) in valid trials");
xlabel("Trials");
ylabel("Mean square error");
ylim([0 25000]);
grid("on");

%%

clc;
for i = 1 % Define training params 

    [trialParams1, trialParams2, trialParams3] = dlDemoThreePattern('xPFC');
    
    outputParams = [{'deepExPFC_V', 1:floor(ModelParametersPFC.NeDeep/3), [200 400] ...
        , 'afr'}; {'deepExPFC_V',ceil(ModelParametersPFC.NeDeep/3):floor(2*ModelParametersPFC.NeDeep/3), ...
        [200 400], 'afr'}; {'deepExPFC_V', ceil(2*ModelParametersPFC.NeDeep/3):ModelParametersPFC.NeDeep, [200 400], 'afr'}; ...
        {'supExPFC_V', 1:ModelParametersPFC.NeSuperficial, [50 700], 'afr'}; ...
        {'midExPFC_V', 1:ModelParametersPFC.NeMid, [50 700], 'afr'}; ...
        {'deepExPFC_V', 1:ModelParametersPFC.NeDeep, [50 700], 'afr'}; ...
        {'supIPVxPFC_V', 1:ModelParametersPFC.NPvSuperficial, [50 700], 'afr'}; ...
        {'midIPVxPFC_V', 1:ModelParametersPFC.NPvMid, [50 700], 'afr'}; ...
        {'deepIPVxPFC_V', 1:ModelParametersPFC.NPvDeep, [50 700], 'afr'}];
    
    targetParams1 = [{'EPenalty', 4:9, 100, 0.01}; {'Compare', [1, 2], 0, 5.5}; {'Compare', [1, 3], 0, 5.5}]; % A 
    targetParams2 = [{'EPenalty', 4:9, 100, 0.01}; {'Compare', [2, 1], 0, 5.5}; {'Compare', [2, 3], 0, 5.5}]; % B
    targetParams3 = [{'EPenalty', 4:9, 100, 0.01}; {'Compare', [3, 1], 0, 5.5}; {'Compare', [3, 2], 0, 5.5}]; % C
    
    dlInputParameters = {trialParams1, trialParams2, trialParams3};
    dlTargetParameters = {targetParams1, targetParams2, targetParams3};
    dlOutputParameters = outputParams;
    
    TBdata = dlTrialBlockGenerator(dlInputParameters, dlTargetParameters, 100, 100);
    
    dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
    dlTrainOptions('dlEpochs') = 10; % % Number of epochs (A.K.A total iterations)
    dlTrainOptions('dlBatchs') = 3; % If a scenario requires the training to be based on a group parameter (e.g mean of errors) use a dlBatch > 1 and set update mode later to batch. 
    dlTrainOptions('dlLambda') = 1e-4; % Higher lambda means more changes based on error, lower may cause model to learn slower or nothing.
        
    dlTrainOptions('dlCheckpoint') = 'true'; % If current step's error is higher based on a threshold, reload last optimal state and continue from that point
    dlTrainOptions('dlCheckpointCoefficient') = 2.047; % A.K.A exploration rate
    dlTrainOptions('dlCheckpointLengthCap') = 11; % If more than 7 steps with no progress passed, return to last checkpoint.
    dlTrainOptions('dlUpdateMode') = 'batch'; % Update on each trial's result or based on batch group results
    
    dlTrainOptions('dlLearningRule') = 'EnhancedDeltaRule'; % Delta rule with a basic change based on biophysical properties 
    dlTrainOptions('dlSimulationFlag') = 1; % If 0, will not run simulations (only for debugging purposes)
    dlTrainOptions('dlOutputLogFlag') = 0; % If 0, will not keep outputs
    dlTrainOptions('dlOfflineOutputGenerator') = 0; % If 1, will generate fake-random outputs (only for debugging purposes)
    
    dlTrainOptions('dlAdaptiveLambda') = 1; % Adaptive lambda parameter; recommended for long simulations.
    dlTrainOptions('dlLambdaCap') = 1.1; % Only if Adaptive lambda is active, recommended to set a upper-bound (UB) or ignore to use default UB (0.01).
    dlTrainOptions('dlExcludeDiverge') = 1; % Exclude non-optimals from model log
    dlTrainOptions('dlTrainExcludeList') = {'Stim'}; % Exclude populations from training
    
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
    
    dlTrainOptions('dlLambda') = 1e-6; % 1e-11(1) -> 1e-4 (4)
    dlTrainOptions('dlAdaptiveLambda') = 0; % Adaptive lambda parameter; recommended for long simulations.
    dlTrainOptions('dlUpdateMode') = 'trial';
    dlTrainOptions('dlLearningRule') = 'BioDeltaRule';
    
    dlTrainOptions('dlTrainExcludeList') = {'Stimuli'};
    dlTrainOptions('dlCheckpointLengthCap') = 20;
    dlTrainOptions('dlEpochs') = 1;
    dlTrainOptions('dlBatchs') = 100;
    
    dlTrainOptions('dlEnhancedMomentum') = 0.4;
    CheckCoeff = 1.28;
    m.dlResetTraining();
    
    dlTrainOptions('dlCustomLog') = ["dlEPowerSpectrum", "dlEPowerSpectrum", "dlEPowerSpectrum", "dlEPowerSpectrum", "dlLFPaxLog", "dlAccuracyBastos2020Task"]; % Name of a function which is in the path
    dlTrainOptions('dlCustomLogArgs') = [argsPowSpectRatio1, argsPowSpectRatio2, argsPowSpectRatio3, argsPowSpectRatio4, argsPowSpectRatio1, argsPowSpectRatio1]; % Arguments of your custom function

end

%%

for  cnt = 1:1

    disp("----------A-----------");
    dlTrainOptions('dlExcludeDiverge') = 1;
    dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff;
    m.dlTrain(TBdata.B1, dlOutputParameters, TBdata.T1, dlTrainOptions);
    
    disp("----------U1----------");
    m.dlErrorsLog = [m.dlErrorsLog, -1];  
    m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 0;
    dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff*1.2;
    m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);
    
    disp("----------B-----------");
    m.dlErrorsLog = [m.dlErrorsLog, -1]; 
    m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 1;
    dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff; 
    m.dlTrain(TBdata.B2, dlOutputParameters, TBdata.T2, dlTrainOptions);
    
    disp("----------U2----------");
    m.dlErrorsLog = [m.dlErrorsLog, -1]; 
    m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 0;
    dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff*1.2; 
    m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);
    
    disp("----------C----------");
    m.dlErrorsLog = [m.dlErrorsLog, -1]; 
    m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 1;
    dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff;
    m.dlTrain(TBdata.B3, dlOutputParameters, TBdata.T3, dlTrainOptions);
    
    disp("----------U3----------");
    m.dlErrorsLog = [m.dlErrorsLog, -1]; 
    m.dlOptimalError = 1e9;dlTrainOptions('dlExcludeDiverge') = 0;
    dlTrainOptions('dlCheckpointCoefficient') = CheckCoeff*1.2;
    m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);

end

m.dlSave();

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
    
    for l = 0:5

        fill([l*w, l*w, l*w+w, l*w+w], [0, errorcap*1.2, errorcap*1.2, 0], [sin(l*0.1), 1, cos(l*0.1)]);hold('on');
        text(l*w+10, errorcap*1.1, tlab(l+1));
        ylim([0 errorcap*1.2]);

    end

end

xlim([0 300]);
plot(x(x>0));grid("on");
title("Errors in trials");

%% Plot Local-field potentials

clc;
m.dlPlotAllPotentials('lfp');

%% Raw run of a simulation (without training)

clc;
tic;
for i = 1:1

    m.dlRunSimulation(dlInputParameters{i}, dlOutputParameters);

end
toc;

%% Beta to Gamma ratio (relative) in trials during predictive task training

clc;

w = 50;
cnt = 1;
tlab = ["A", "U1", "B", "U2", "C", "U3", "A", "U4", "B", "U5", "C", "U6"];
f = figure("Position", [0 0 1500 1000]);

set1 = [2:5, 16:19];
set2 = [6:9, 20:23];

ratioLog = {m.dlCustomLog(1, :)};
ratioLog = cell2mat(ratioLog{1});

for i = 1:4

    for l = 0:5

        subplot(2, 2, i);
        fill([l*w, l*w, l*w+w, l*w+w], [0, 40, 40, 0], [sin(l*0.1), 1, cos(l*0.1)], 'HandleVisibility','off');
        hold('on');

    end

end

dispLabels = m.dlCustomLogLabel{1};

for i = 2:9

    subplot(2, 2, mod(cnt-1, 4)+1);
    plot(ratioLog(i, :), 'DisplayName', dispLabels{i});
    legend("Location", "southwest");

    if i > 4

        q = ratioLog(i-4, :);
        q2 = ratioLog(i, :);
        ql = min(min(q), min(q2));
        qu = max(max(q), max(q2));
    
        ylim([ql-1 qu+1]);
        title(dispLabels{i} + " and " + dispLabels{i-4});

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