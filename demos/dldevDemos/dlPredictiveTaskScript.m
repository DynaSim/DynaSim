%% Simulation of predictive task implemented on DynaSim/DynaLearn via Reinforcement learning 

% Based on [Bastos2020:Layer and rhythm specificity for predictive routing, A.M Bastos et. al, 2020]
% Requirement: DynaSim/dev latest version
% Minimum MatLab version: R2019a
% 5th of July, 2022
% @HNXJ

%% Model parameters

clear;clc;

%%% Area PFC layer sizes (relative)
NeS_PFC = 71;NiS_PFC = 14;
NeM_PFC = 47;NiM_PFC = 7;
NeD_PFC = 47;NiD_PFC = 4;

Nin_PFC = 6;
Nout_PFC = 6;
NoiseRate_PFC = 7;
Nstim_PFC = 3;

%%% Area V4 layer sizes (relative) 
NeS_V4 = 77;NiS_V4 = 17;
NeM_V4 = 49;NiM_V4 = 11;
NeD_V4 = 41;NiD_V4 = 7;

Nin_V4 = 6;
Nout_V4 = 6;
NoiseRate_V4 = 4;
Nstim_V4 = 3;

%%% Call Laminar Cortex Constructor Functions
dsCellPFC = dlLaminarCortexNet(NeS_PFC, NiS_PFC, NeM_PFC, NiM_PFC, NeD_PFC, ...
    NiD_PFC, Nin_PFC, Nout_PFC, Nstim_PFC, NoiseRate_PFC, 'PFC'); % Laminar PFC model with specific parameters
dsCellV4 = dlLaminarCortexNet(NeS_V4, NiS_V4, NeM_V4, NiM_V4, NeD_V4, ...
    NiD_V4, Nin_V4, Nout_V4, Nstim_V4, NoiseRate_V4, 'V4'); % Laminar V4 model with specific parameters

%%% Connect two models; consider adding connections later (TODO:ADD_CONN)
dsModel = dsCombineSpecifications(dsCellV4, dsCellPFC); % Full model

%%% Connecting two models
connectionWeigth1 = 0.21*rand(NeD_V4, NeM_PFC) + 0.27;
c = length(dsModel.connections) + 1;
connection1.direction = [dsCellV4.populations(5).name, '->', dsCellPFC.populations(3).name];

connection1.source = dsCellV4.populations(5).name;
connection1.target = dsCellPFC.populations(3).name;
connection1.mechanism_list={'iPoisson'};
connection1.parameters={'gAMPA', .24, 'tauAMPA', 4.17, 'netcon', connectionWeigth1};

connectionWeigth2 = 0.27*rand(NiD_V4, NiM_PFC) + 0.21;
c = length(dsModel.connections) + 1;
connection2.direction = [dsCellV4.populations(6).name, '->', dsCellPFC.populations(4).name];

connection2.source = dsCellV4.populations(6).name;
connection2.target = dsCellPFC.populations(4).name;
connection2.mechanism_list={'iPoisson'};
connection2.parameters={'gAMPA', .27, 'tauAMPA', 3.74, 'netcon', connectionWeigth2};

dsModel = dlConnectModels({dsCellV4, dsCellPFC}, {connection1, connection2});

%% Create DynaLearn Class (Only first time, if file does not exist already)

tic;m = DynaLearn(dsModel, 'models/dlModelPredictivePFC', 'mex');toc; % ~70 min, MEXGEN or ~1 min, RAWGEN
m.dlSave(); % < 1sec

%% Load DynaLearn Class

clear;clc;

m = DynaLearn(); % ~ 1sec
% m = m.dlLoad('models/dlModelPredictivePFC3'); % ~ 10sec, Trained for ~1200 trials
% m = m.dlLoad('models/dlModelPredictivePFC4'); % ~ 10sec, New! keeping track of its activity in Gamma/Beta **
m = m.dlLoad('models/dlModelPredictivePFC5'); % ~ 10sec, New larger model; keeping track of its activity in Gamma/Beta **

% m.dlSimulate(); % ~ 40sec

%% Trial: training script preparation, 50-block and 50-trial

[trialParams1, trialParams2, trialParams3] = dlDemoThreePattern();

outputParams = [{'DeepE_V', 1:4, [300 500], 'afr'}; {'DeepE_V', 5:8, [300 500], 'afr'}; {'DeepE_V', 9:12, [300 500], 'afr'}];
targetParams1 = [{'MSE', 1, 27, 0.2}; {'MSE', 2, 20, 0.2}; {'MSE', 3, 20, 0.2}; {'Compare', [1, 2], 0, 0.3}; {'Compare', [1, 3], 0, 0.3}; {'Diff', [2, 3], 0, 0.04}]; % A 
targetParams2 = [{'MSE', 2, 27, 0.2}; {'MSE', 1, 20, 0.2}; {'MSE', 3, 20, 0.2}; {'Compare', [2, 1], 0, 0.3}; {'Compare', [2, 3], 0, 0.3}; {'Diff', [1, 3], 0, 0.04}]; % B
targetParams3 = [{'MSE', 3, 27, 0.2}; {'MSE', 2, 20, 0.2}; {'MSE', 1, 20, 0.2}; {'Compare', [3, 1], 0, 0.3}; {'Compare', [3, 2], 0, 0.3}; {'Diff', [1, 2], 0, 0.04}]; % C

dlInputParameters = {trialParams1, trialParams2, trialParams3};
dlTargetParameters = {targetParams1, targetParams2, targetParams3};
dlOutputParameters = outputParams;

TBdata = dlTrialBlockGenerator(dlInputParameters, dlTargetParameters, 200, 2000);

dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
dlTrainOptions('dlEpochs') = 500; % % Number of epochs (A.K.A total iterations)
dlTrainOptions('dlBatchs') = 3; % If a scenario requires the training to be based on a group parameter (e.g mean of errors) use a dlBatch > 1 and set update mode later to batch. 
dlTrainOptions('dlLambda') = 1e-5; % Higher lambda means more changes based on error, lower may cause model to learn slower or nothing.
    
dlTrainOptions('dlCheckpoint') = 'true'; % If current step's error is higher based on a threshold, reload last optimal state and continue from that point
dlTrainOptions('dlCheckpointCoefficient') = 1.74; % A.K.A exploration rate
dlTrainOptions('dlBadTrialEliminatorFlag') = 1; % A.K.A backtrack of only usefull trials
dlTrainOptions('dlUpdateMode') = 'batch'; % Update on each trial's result or based on batch group results

dlTrainOptions('dlLearningRule') = 'BioDeltaRule'; % Delta rule with a basic change based on biophysical properties 
dlTrainOptions('dlSimulationFlag') = 1; % If 0, will not run simulations (only for debugging purposes)
dlTrainOptions('dlOutputLogFlag') = 1; % If 0, will not keep outputs
dlTrainOptions('dlOfflineOutputGenerator') = 0; % If 1, will generate fake-random outputs (only for debugging purposes)

dlTrainOptions('dlAdaptiveLambda') = 1; % Adaptive lambda parameter; recommended for long simulations.
dlTrainOptions('dlLambdaCap') = 3e-2; % Only if Adaptive lambda is active, recommended to set a upper-bound (UB) or ignore to use default UB (0.01).

% Initial training on the model to reach a plausible local minimia like
% the task in the paper the model should also learn the basics of the task.
% We shortly train the model by cues to put it close to a local minimia.

dlTrainOptions('dlLambda') = 7e-5;
dlTrainOptions('dlEpochs') = 3;
dlTrainOptions('dlBatchs') = 3;

argsPSR = struct();

argsPSR.lf1 = 7;
argsPSR.hf1 = 31;
argsPSR.lf2 = 47;
argsPSR.hf2 = 74;

dlTrainOptions('dlCustomLog') = "dlPowerSpectrumRatio"; % Name of a function which is in the path
dlTrainOptions('dlCustomLogArgs') = argsPSR; % Arguments of your custom function

%%

clc;
m.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

%%

clc;
figure();
tlab = ["A", "Unpred.", "B", "Unpred", "C", "Unpred"];
w = 200;

for k = 1:3
    
    subplot(3, 1, k);
    for l = 0:5
        fill([l*w, l*w, l*w+w, l*w+w], [0, 20, 20, 0], [1, cos(l*0.15), sin(l*0.15)], 'DisplayName', tlab(l+1));hold('on');
        text(l*w+10, 19, tlab(l+1));
    end
    
    plot(m.dlCustomLog(k*2-1, :), 'DisplayName', 'E');
    plot(m.dlCustomLog(k*2, :), 'DisplayName', 'I');
    xlim([0, w*6.5]);
    xlabel("$S(\frac{\beta}{\gamma}) \: of \: " + m.dlCustomLogLabel(k*2) + "\: \& \:" + m.dlCustomLogLabel(k*2+1) + "$", 'Interpreter', 'latex');
    grid("on");legend();
end

%% Block-trial phase

clc;
dlTrainOptions('dlLambda') = 1e-4;
dlTrainOptions('dlUpdateMode') = 'trial';
dlTrainOptions('dlEpochs') = 1;
dlTrainOptions('dlBatchs') = 200;

argsPSR = struct();

argsPSR.lf1 = 8;
argsPSR.hf1 = 32;
argsPSR.lf2 = 40;
argsPSR.hf2 = 100;

dlTrainOptions('dlCustomLog') = "dlPowerSpectrumRatio"; % Name of a function which is in the path
dlTrainOptions('dlCustomLogArgs') = argsPSR; % Arguments of your custom function

dlTrainOptions('dlCheckpointCoefficient') = 2.4;
m.dlTrain(TBdata.B1, dlOutputParameters, TBdata.T1, dlTrainOptions);

m.dlOptimalError = 1e9;
dlTrainOptions('dlCheckpointCoefficient') = 4.7;
m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);

m.dlOptimalError = 1e9;
dlTrainOptions('dlCheckpointCoefficient') = 2.4; 
m.dlTrain(TBdata.B2, dlOutputParameters, TBdata.T2, dlTrainOptions);

m.dlOptimalError = 1e9;
dlTrainOptions('dlCheckpointCoefficient') = 4.7; 
m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);

m.dlOptimalError = 1e9;
dlTrainOptions('dlCheckpointCoefficient') = 2.4;
m.dlTrain(TBdata.B3, dlOutputParameters, TBdata.T3, dlTrainOptions);

m.dlOptimalError = 1e9;
dlTrainOptions('dlCheckpointCoefficient') = 4.7;
m.dlTrain(TBdata.TrB, dlOutputParameters, TBdata.TrT, dlTrainOptions);

%% Errors log plot

% clc;

wk = 4;
figure('position', [0, 0, 1400, 700]);
n = max(size(m.dlErrorsLog));
x = zeros(1, ceil(n/wk));

for i = 0:wk:n-wk
    x(ceil((i+1)/wk)) = mean(m.dlErrorsLog(i+1:i+wk));
end

w = 50;
for k = 1:6
    
    for l = 0:5
        fill([l*w, l*w, l*w+w, l*w+w], [0, 100, 100, 0], [sin(l*0.2), 1, cos(l*0.2)]);hold('on');
        text(l*w+10, 89, tlab(l+1));
    end
end

plot(x);grid("on");
title("Errors in batchs");
% m.dlPlotBatchErrors(2);

%% Plot Local-field potentials

clc;
m.dlPlotAllPotentials('lfp');

%% Run a simulation (without training)

for i = 1:3
    m.dlRunSimulation(dlInputParameters{i}, dlOutputParameters);
%     m.dlPlotAllPotentials('lfp');
end

%%

clc;
opts = containers.Map();
% opts("lf") = 50;
% opts("hf") = 100;
% m.dlPlotAllPotentials('avgfft', opts);

opts("lf") = 50;
opts("hf") = 100;
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

%% G/B ratio log

dtf = ceil(1 / (obj.dldT*obj.dlDownSampleFactor));

lf = opts("lf")*dtf;
hf = opts("hf")*dtf; 
freqCap = 0;

for i = (k-1)*6+1:min((k*6), 6)

    x = dlPotentials{1, i+1};
    fqs = linspace(1, 500, max(size(x)));
    subplot((min(k*6, n-1) - (k-1)*6), 1, mod(i-1, (min(k*6, n-1) - (k-1)*6))+1);
    ffts = abs(fft(mean(x, 2))) * min(size(x)) / 1000;
    yf = smooth(ffts(lf:hf));
    area(fqs(lf:hf), yf);grid("on");

    if freqCap == 0
        freqCap = max(ffts(lf:hf))*1.2;
        ylim([0, freqCap]);
    else
        ylim([0, freqCap]);
    end

    ylabel(dlLabels(i+1));

end

disp("Temp edit for 6 subplots; average fft");
xlabel(mode + " in frequency (Hz)");
                    
%% End