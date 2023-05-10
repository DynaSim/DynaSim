%% Simulation of predictive task implemented on DynaSim/DynaLearn via Reinforcement learning 

% Requirement: DynaSim/dev_2022 latest version
% Minimum MatLab version: R2021a
% @April2023
% @HNXJ

%% Initiation

clear;clc;
cd('D:\Works\Computational');
addpath(genpath('DynaSim'));
cd('DynaSim');

%% AutoRunSc: Pre-fitting

clear;clc;

CurrentSize = 64;
TotalSize = ones(1, 20)*CurrentSize;
noise_rate = 7.4;
performance_coefficient = 0; % Pre-fit

ResetOptimalError = 'on';
RemakeFlag = 1;
tune_flag = 1;
epochs = 100; % Low number for implementation and debugging purposes

%%% Create model parameters struct example
% ModelName = "V1";
% ModelName = "MST";
ModelName = "PFC";
ModelParameters = struct();

%%% Area PFC layer sizes (relative)
ModelParameters.NeSuperficial = ceil(0.21*CurrentSize);
ModelParameters.NPvSuperficial = ceil(0.02*CurrentSize);
ModelParameters.NSomSuperficial = ceil(0.04*CurrentSize);
ModelParameters.NVipSuperficial = ceil(0.03*CurrentSize);

ModelParameters.NeMid = ceil(0.15*CurrentSize);
ModelParameters.NPvMid = ceil(0.03*CurrentSize);
ModelParameters.NSomMid = ceil(0.01*CurrentSize);
ModelParameters.NVipMid = ceil(0.01*CurrentSize);

ModelParameters.NeDeep = ceil(0.41*CurrentSize);
ModelParameters.NPvDeep = ceil(0.03*CurrentSize);
ModelParameters.NSomDeep = ceil(0.03*CurrentSize);
ModelParameters.NVipDeep = ceil(0.03*CurrentSize);

ModelParameters.Nin = 6;
ModelParameters.Nout = 6;
ModelParameters.NoiseRate = noise_rate; % 10%
ModelParameters.Nstim = 3;

% dsCellLaminar = dlLaminarCortexNetNL(ModelParameters, ModelName);

for model_size_id = 1:1

    CurrentSize = TotalSize(model_size_id);
    dlPassiveDynamicsPerformer(RemakeFlag, ResetOptimalError, ModelName, ModelParameters, model_size_id, performance_coefficient, tune_flag, epochs);

end

%% AutoRunScript: Main simulation

clear;clc;

TotalSize = ones(1, 20)*64;
noise_rate = 7.4;
performance_coefficient = 1; % Main task

ResetOptimalError = 'on';
RemakeFlag = 0;
tune_flag = 0;

for model_size_id = 1:1

    CurrentSize = TotalSize(model_size_id);
    dlPassiveDynamicsPerformer(RemakeFlag, ResetOptimalError, ModelName, ModelParameters, model_size_id, performance_coefficient, tune_flag);

end

%% Plot sample response (loader)

clear;clc;
id = 1;
m = DynaLearn(); % ~ 1sec
m = m.dlLoad(char("dlModels/dlPredictiveCorticalCircuitModelNL" + string(id))); % ~ 10sec, New larger model; keeping track of its activity in Gamma/Beta **

%% Plots
tic;
trialParamsTemp = containers.Map();
trialParamsTemp('tspan') = [0 10000];
m.dlUpdateParams(trialParamsTemp);
m.dlSimulate();

opts = containers.Map();
opts("lf") = 1;
opts("hf") = 60;
m.dlPlotAllPotentials('raster');
toc;
% m.dlPlotAllPotentials('avgfft', opts);
% opts = struct();opts.name = "test1";
% m.dlPlotAllPotentials('lfp', opts);  
% m.dlPlotAllPotentials('lfpsave', opts);   
% m.dlPlotAllPotentials('lfpmap', opts);   

%% TB1

%%

figure('Position', [0, 0, 1700, 1400]);

for i = 5

    load("dlModels/dlTCTModels2/solve/signals_test" + num2str(i) + ".mat");

%     subplot(3, 3, i);
    imagesc(imgRaster>4);
    colormap("gray");

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