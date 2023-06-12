% addpath(genpath('/Users/jason/Documents/me/docs - research/andre/BetaGammaPushPull/DynaSim'));

clear;close all;clc;
PathToDynaSim = 'D:\Works\Computational'; % Change it based on your local path for dynasim
cd(PathToDynaSim);
addpath(genpath('DynaSim'));
cd('DynaSim');

%% 1. PING/PINB only

% 1. Optimize '(INfast->INfast, INfast->ES)','tauD' to achieve fnat=15,20,25Hz
%   This will become our tauGABAslow
% 2. Optimize '(INfast->INfast, INfast->ES)','tauD' to achieve fnat=40,50,60Hz
%   This will become our tauGABAfast

% define equations of cell model (same for E and I populations)

eqns={
  'dV/dt=Iapp+@current+noise*randn(1,N_pop); Iapp=0; noise=0'
};

Iapp = 10;

gGABAfast = .1;
gAMPA = .1;

tauGABAfast = 5;
tauAMPA = 2;

s=[];
s.populations(1).name='ES';
s.populations(1).size=80;
s.populations(1).equations=eqns;
s.populations(1).mechanism_list={'iNa','iK'};
s.populations(1).parameters={'Iapp',Iapp,'gNa',120,'gK',36,'noise', 40};

s.populations(2).name='INfast';
s.populations(2).size=20;
s.populations(2).equations=eqns;
s.populations(2).mechanism_list={'iNa','iK'};
s.populations(2).parameters={'Iapp',0,'gNa',120,'gK',36,'noise', 40};

s.connections(1).direction='INfast->ES';
s.connections(1).mechanism_list={'iGABAa'};
s.connections(1).parameters={'tauD',tauGABAfast,'gGABAa',gGABAfast};

s.connections(2).direction='ES->INfast';
s.connections(2).mechanism_list={'iAMPA'};
s.connections(2).parameters={'tauD',tauAMPA,'gAMPA',gAMPA};

s.connections(3).direction='INfast->INfast';
s.connections(3).mechanism_list={'iGABAa'};
s.connections(3).parameters={'tauD',tauGABAfast,'gGABAa',gGABAfast};

% vary = {'(INfast->INfast, INfast->ES)','tauD',[1 5 10 20 50]};
% 
% tic; data = dsSimulate(s, 'tspan', [0 1000], 'downsample_factor', 10, ...
%                      'solver', 'euler', 'dt', .002, ...
%                      'compile_flag', 0, 'parallel_flag', 1, ...
%                      'vary', vary, 'verbose_flag', 1); toc;
% 
% dsPlot(data);
% dsPlot(data, 'plot_type','raster');
% tic; dsPlot(data, 'plot_type','power'); toc;

%% 1.0 Optimize firing rate

RunDuration = 500;
study_dir = 'dlModels/EI4';
target_FRQ = 70; % Hz
target_AV = -58; % mV
dl = DynaLearn(s, study_dir, 'mex', 'EI', 1);
dl.dlSave();

%%

dl.dlLoader(study_dir);

%% Define optimization parameters

dlInputParameters = {dlNullInputs(RunDuration)}; % Null inputs (no external stimulation).
dlOutputParameters = [{['ES', '_V'], 1:80, [10 RunDuration], 'afr'}; {['ES', '_V'], 1:80, [10 RunDuration], 'av'}];
dlTargetParameters = {[{'MSE', 1, target_FRQ, 1}; {'MSE', 2, -58, 1}]}; % Format: (1)Mode;(2)Output indices;(3)Target value;(4)Weight in loss equation

% Define training options
dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
dlTrainOptions('dlLambda') = 1e-5; % Fitting/Learning rate
dlTrainOptions('dlEpochs') = 300; % Total iterations 
dlTrainOptions('dlEnhancedMomentum') = 0.4;

dlTrainOptions('dlExcludeDiverge') = 0; % Exclude non-optimals from model log
dlTrainOptions('dlCheckpointLengthCap') = 11;
dlTrainOptions('dlUpdateVerbose') = 1;
dlTrainOptions('dlTrainIncludeList') = ["ES_INfast_iGABAa_tauD", "ES_noise"];

dlTrainOptions('dlTrainRestrictList') = ["ES_noise", "INfast_noise"]; % Restrict list
dlTrainOptions('dlTrainRestrictCoef') = {.001, .001}; % Restrict coeffs
dlTrainOptions('dlTrainSyncList') = [["ES_INfast_iGABAa_tauD", ...
    "INfast_INfast_iGABAa_tauD"]; ["ES_noise", "INfast_noise"]]; % Each array row will be overrided
                                  % by synchronizing to the first element
dlTrainOptions('dlCheckpointCoefficient') = 1.5;

dl.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

%%

dl.dlLoadOptimal 
dl.dlSimulate
dl.dlPlotAllPotentials('lfp');
dl.dlPlotBatchErrors();

%%

p = load([study_dir, '/Optimalparams.mat']);
p = p.p;
fprintf("Starting parameters > ES_noise = %f, INfast_noise = %f \n ->ES_INfast_iGABAa_tauD = %f\n->INfast_INfast_iGABAa_tauD = %f\n", 40.0, 40.0, 5.0, 5.0);
fprintf("Optimal parameters > ES_noise = %f, INfast_noise = %f \n ->ES_INfast_iGABAa_tauD = %f\n->INfast_INfast_iGABAa_tauD = %f\n", p.ES_noise, p.INfast_noise, p.ES_INfast_iGABAa_tauD, p.INfast_INfast_iGABAa_tauD);
fprintf("Optimal outputs > ES_V_aFR = %f, ES_V_aV = %f\n", dl.dlLastOutputs{1}, dl.dlLastOutputs{2});
fprintf("Target outputs > ES_V_aFR = %f, ES_V_aV = %f\n", target_FRQ, target_AV);

%% 1.1 Optimize natural frequency (Gamma~50Hz)

RunDuration = 500;
study_dir = 'dlModels/EI5';
target_FRQ = 50; % Hz
target_AV = -58; % mV
dl = DynaLearn(s, study_dir, 'mex', 'EI', 1);
dl.dlSave();

%%

dl.dlLoader(study_dir);

%% Define optimization parameters

dlInputParameters = {dlNullInputs(RunDuration)}; % Null inputs (no external stimulation).
dlOutputParameters = [{['ES', '_V'], 1:80, [10 RunDuration], 'fnat'}; {['ES', '_V'], 1:80, [10 RunDuration], 'av'}; {['ES', '_V'], 1:80, [10 RunDuration], 'afr'}];
dlTargetParameters = {[{'MSE', 1, target_FRQ, 1}; {'MSE', 2, -58, 1}]}; % Format: (1)Mode;(2)Output indices;(3)Target value;(4)Weight in loss equation

% Define training options
dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
dlTrainOptions('dlLambda') = 1e-3; % Fitting/Learning rate
dlTrainOptions('dlEpochs') = 1; % Total iterations 
dlTrainOptions('dlEnhancedMomentum') = 0.2;

dlTrainOptions('dlExcludeDiverge') = 0; % Exclude non-optimals from model log
dlTrainOptions('dlCheckpointLengthCap') = 11;
dlTrainOptions('dlUpdateVerbose') = 1;
dlTrainOptions('dlTrainIncludeList') = ["ES_INfast_iGABAa_tauD", "ES_noise"];

dlTrainOptions('dlTrainRestrictList') = ["ES_noise", "INfast_noise"]; % Restrict list
dlTrainOptions('dlTrainRestrictCoef') = {.001, .001}; % Restrict coeffs
dlTrainOptions('dlTrainSyncList') = [["ES_INfast_iGABAa_tauD", ...
    "INfast_INfast_iGABAa_tauD"]; ["ES_noise", "INfast_noise"]]; % Each array row will be overrided
                                  % by synchronizing to the first element
dlTrainOptions('dlCheckpointCoefficient') = 1.5;

dl.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

%%

dl.dlLoadOptimal 
dl.dlSimulate
dl.dlPlotAllPotentials('lfp');
dl.dlPlotBatchErrors();

%%

p = load([study_dir, '/Optimalparams.mat']);
p = p.p;
fprintf("Starting parameters > ES_noise = %f, INfast_noise = %f \n ->ES_INfast_iGABAa_tauD = %f\n->INfast_INfast_iGABAa_tauD = %f\n", 40.0, 40.0, 5.0, 5.0);
fprintf("Optimal parameters > ES_noise = %f, INfast_noise = %f \n ->ES_INfast_iGABAa_tauD = %f\n->INfast_INfast_iGABAa_tauD = %f\n", p.ES_noise, p.INfast_noise, p.ES_INfast_iGABAa_tauD, p.INfast_INfast_iGABAa_tauD);
fprintf("Optimal outputs > ES_V_aFR = %f, ES_V_aV = %f\n", dl.dlLastOutputs{1}, dl.dlLastOutputs{2});
fprintf("Target outputs > ES_V_aFR = %f, ES_V_aV = %f\n", target_FRQ, target_AV);

%% 1.2 Optimize natural frequency (Beta~20Hz)

RunDuration = 500;
study_dir = 'dlModels/EI6';
target_FRQ = 20; % Hz
target_AV = -70; % mV
dl = DynaLearn(s, study_dir, 'mex', 'EI', 1);
dl.dlSave();

%%

dl.dlLoader(study_dir);

%% Define optimization parameters

dlInputParameters = {dlNullInputs(RunDuration)}; % Null inputs (no external stimulation).
dlOutputParameters = [{['ES', '_V'], 1:80, [10 RunDuration], 'fnat'}; {['ES', '_V'], 1:80, [10 RunDuration], 'av'}; {['ES', '_V'], 1:80, [10 RunDuration], 'afr'}];
dlTargetParameters = {[{'MQE', 1, target_FRQ, 1}; {'MSE', 2, target_AV, 1}]}; % Format: (1)Mode;(2)Output indices;(3)Target value;(4)Weight in loss equation

% Define training options
dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
dlTrainOptions('dlLambda') = 1e-2; % Fitting/Learning rate
dlTrainOptions('dlEpochs') = 400; % Total iterations 
dlTrainOptions('dlEnhancedMomentum') = 0.5;

dlTrainOptions('dlExcludeDiverge') = 0; % Exclude non-optimals from model log
dlTrainOptions('dlCheckpointLengthCap') = 11;
dlTrainOptions('dlUpdateVerbose') = 1;
dlTrainOptions('dlTrainIncludeList') = ["ES_INfast_iGABAa_tauD", "ES_noise"];

dlTrainOptions('dlTrainRestrictList') = ["ES_noise", "INfast_noise"]; % Restrict list
dlTrainOptions('dlTrainRestrictCoef') = {.001, .001}; % Restrict coeffs
dlTrainOptions('dlTrainSyncList') = [["ES_INfast_iGABAa_tauD", ...
    "INfast_INfast_iGABAa_tauD"]; ["ES_noise", "INfast_noise"]]; % Each array row will be overrided
                                  % by synchronizing to the first element
dlTrainOptions('dlCheckpointCoefficient') = 1.5;

dl.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

%%

dl.dlLoadOptimal 
dl.dlSimulate
dl.dlPlotAllPotentials('lfp');
dl.dlPlotBatchErrors();

%%

p = load([study_dir, '/Optimalparams.mat']);
p = p.p;
fprintf("Starting parameters > ES_noise = %f, INfast_noise = %f \n ->ES_INfast_iGABAa_tauD = %f\n->INfast_INfast_iGABAa_tauD = %f\n", 40.0, 40.0, 5.0, 5.0);
fprintf("Optimal parameters > ES_noise = %f, INfast_noise = %f \n ->ES_INfast_iGABAa_tauD = %f\n->INfast_INfast_iGABAa_tauD = %f\n", p.ES_noise, p.INfast_noise, p.ES_INfast_iGABAa_tauD, p.INfast_INfast_iGABAa_tauD);
fprintf("Optimal outputs > ES_V_fNAT = %f, ES_V_aV = %f\n", dl.dlLastOutputs{1}, dl.dlLastOutputs{2});
fprintf("Target outputs > ES_V_fNAT = %f, ES_V_aV = %f\n", target_FRQ, target_AV);

%% 1.3 Optimize natural frequency (Low-gamma~40Hz)

RunDuration = 500;
study_dir = 'dlModels/EI7';
target_FRQ = 20; % Hz
target_AV = -70; % mV
dl = DynaLearn(s, study_dir, 'mex', 'EI', 1);
dl.dlSave();

%%

dl = DynaLearn.dlLoader(study_dir);

%% Define optimization parameters

dlInputParameters = {dlNullInputs(RunDuration)}; % Null inputs (no external stimulation).
dlOutputParameters = {{['ES', '_V'], 1:80, [10 RunDuration], 'fnat'}; {['ES', '_V'], 1:80, [10 RunDuration], 'av'}; {['INfast', '_V'], 1:20, [10 RunDuration], 'fnat'}};
dlTargetParameters = {{'RPenalty', 1:2, 10, 1, 8, 30, 40, 90}; {'MSE', 1, target_FRQ, 1}; {'MAE', 2, -65, 1}; {'MSE', 3, target_FRQ, 1}}; % Format: (1)Mode;(2)Output indices;(3)Target value;(4)Weight in loss equation

% Define training options
dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
dlTrainOptions('dlLambda') = 1e-3;% Fitting/Learning rate
dlTrainOptions('dlEpochs') = 2; % Total iterations 
dlTrainOptions('dlEnhancedMomentum') = 0.25;

dlTrainOptions('dlExcludeDiverge') = 0; % Exclude non-optimals from model log
dlTrainOptions('dlCheckpointLengthCap') = 11;
dlTrainOptions('dlUpdateVerbose') = 1;
dlTrainOptions('dlTrainIncludeList') = ["ES_INfast_iGABAa_tauD", "ES_noise", "INfast_INfast_iGABAa_tauD"];

dlTrainOptions('dlTrainRestrictList') = ["ES_noise", "INfast_noise"]; % Restrict list
dlTrainOptions('dlTrainRestrictCoef') = {.001, .001}; % Restrict coeffs
dlTrainOptions('dlTrainSyncList') = [["ES_INfast_iGABAa_tauD", "INfast_INfast_iGABAa_tauD"]; ["ES_noise", "INfast_noise"]]; % Each array row will be overrided
                                  % by synchronizing to the first element
dlTrainOptions('dlCheckpointCoefficient') = 4.5;

dl.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

%%

dl.dlLoadOptimal 
dl.dlSimulate
dl.dlPlotAllPotentials('lfp');
dl.dlPlotBatchErrors();

%%

p = load([study_dir, '/Optimalparams.mat']);
p = p.p;
fprintf("Starting parameters > ES_noise = %f, INfast_noise = %f \n ->ES_INfast_iGABAa_tauD = %f\n->INfast_INfast_iGABAa_tauD = %f\n", 40.0, 40.0, 5.0, 5.0);
fprintf("Optimal parameters > ES_noise = %f, INfast_noise = %f \n ->ES_INfast_iGABAa_tauD = %f\n->INfast_INfast_iGABAa_tauD = %f\n", p.ES_noise, p.INfast_noise, p.ES_INfast_iGABAa_tauD, p.INfast_INfast_iGABAa_tauD);
fprintf("Optimal outputs > ES_V_fNAT = %f, ES_V_aV = %f\n", dl.dlLastOutputs{1}, dl.dlLastOutputs{2});
fprintf("Target outputs > ES_V_fNAT = %f, ES_V_aV = %f\n", target_FRQ, target_AV);


%% 2. Optimize contextual input for beta/gamma push-pull
% define equations of cell model (same for E and I populations)

eqns={
  'HH: dV/dt=IappCue+IappContext.*Context+@current+noise*randn(1,N_pop); Context=1; IappCue=0; IappContext=0; noise=0'
};

IappCue = 10;
IappContext = 10;
gGABAslow = .1;
gGABAfast = .1;
gAMPA = .1;

tauGABAslow = 20;
tauGABAfast = 5;
tauAMPA = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=[];
s.populations(1).name='ES';
s.populations(1).size=80;
s.populations(1).equations=eqns;
s.populations(1).mechanism_list={'iNa','iK'};
s.populations(1).parameters={'IappCue',IappCue, 'IappContext', IappContext,'gNa',120,'gK',36,'noise',40};
s.populations(2).name='INfast';
s.populations(2).size=20;
s.populations(2).equations=eqns;
s.populations(2).mechanism_list={'iNa','iK'};
s.populations(2).parameters={'IappCue',IappCue, 'IappContext', IappContext,'gNa',120,'gK',36,'noise',40};
s.populations(3).name='INslow';
s.populations(3).size=20;
s.populations(3).equations=eqns;
s.populations(3).mechanism_list={'iNa','iK'};
s.populations(3).parameters={'IappCue',0,'gNa',120,'gK',36,'noise',40};
s.connections(1).direction='INfast->ES';
s.connections(1).mechanism_list={'iGABAa'};
s.connections(1).parameters={'tauD',tauGABAfast,'gGABAa',gGABAfast};
s.connections(2).direction='ES->INfast';
s.connections(2).mechanism_list={'iAMPA'};
s.connections(2).parameters={'tauD',tauAMPA,'gAMPA',gAMPA};
s.connections(3).direction='INslow->ES';
s.connections(3).mechanism_list={'iGABAa'};
s.connections(3).parameters={'tauD',tauGABAslow,'gGABAa',gGABAslow};
s.connections(4).direction='ES->INslow';
s.connections(4).mechanism_list={'iAMPA'};
s.connections(4).parameters={'tauD',tauAMPA,'gAMPA',gAMPA};
s.connections(5).direction='INslow->INfast';
s.connections(5).mechanism_list={'iGABAa'};
s.connections(5).parameters={'tauD',tauGABAslow,'gGABAa',gGABAslow};
s.connections(6).direction='INslow->INslow';
s.connections(6).mechanism_list={'iGABAa'};
s.connections(6).parameters={'tauD',tauGABAslow,'gGABAa',gGABAslow};
s.connections(7).direction='INfast->INslow';
s.connections(7).mechanism_list={'iGABAa'};
s.connections(7).parameters={'tauD',tauGABAfast,'gGABAa',gGABAfast};
s.connections(8).direction='INfast->INfast';
s.connections(8).mechanism_list={'iGABAa'};
s.connections(8).parameters={'tauD',tauGABAfast,'gGABAa',gGABAfast};

% vary = {'INfast','IappCue',[0 IappCue]; 'INslow','IappContext',[0 IappContext/2 IappContext]};
% tic
% data = dsSimulate(s, 'tspan', [0 1000], 'downsample_factor', 10, ...
%                      'solver', 'euler', 'dt', .002, ...
%                      'compile_flag', 0, 'parallel_flag', 1, ...
%                      'vary', vary, 'verbose_flag', 1);
% toc
% ​
% dsPlot(data);
% dsPlot(data, 'plot_type','raster');
% tic; dsPlot(data, 'plot_type','power'); toc

%%

RunDuration = 500;
study_dir = 'dlModels/EIX3';
target_FRQ = 20; % Hz
target_AV = -70; % mV
dl = DynaLearn(s, study_dir, 'mex', 'EI', 1);
dl.dlSave();

%% If already exists

dl = DynaLearn.dlLoader(study_dir);

%% Define optimization parameters


% Optimize IappContext to ES, INfast, INslow 
% Condition 1 (context=0): gamma>beta
% Condition 2 (context=1): beta>gamma
%   where gamma = 55-75Hz, beta = 20-40Hz
% ​
% to ES, INfast, INslow
% ​
% targets (Rpenalty)
% condition1: gamma > beta
% condition2: beta > gamma

condition1 = containers.Map();
condition2 = containers.Map();
condition1('ES_Context') = 0;
condition1('INfast_Context') = 0;
condition1('INslow_Context') = 0;
condition1('tspan') = [0 RunDuration];
condition2('ES_Context') = 1;
condition2('INfast_Context') = 1;
condition2('INslow_Context') = 1;
condition2('tspan') = [0 RunDuration];

dlInputParameters = {condition1, condition2}; % Context inputs (no external stimulation).
dlOutputParameters = {{['ES', '_V'], 1:80, [10 RunDuration], 'av'}; {['INfast', '_V'], 1:20, [10 RunDuration], 'av'}; {['INslow', '_V'], 1:20, [10 RunDuration], 'av'}};
dlTargetParameters1 = {{'QRPenalty', 1, 5, 1, 20, 40, 55, 75}; {'MAE', 1, -70, 1}; {'MAE', 2, -70, 1}; {'MAE', 3, -70, 1}}; % Format: (1)Mode;(2)Output indices;(3)Target value;(4)Weight in loss equation
dlTargetParameters2 = {{'QRPenalty', 1, .2, 1, 20, 40, 55, 75}; {'MAE', 1, -70, 1}; {'MAE', 2, -70, 1}; {'MAE', 3, -70, 1}}; % Format: (1)Mode;(2)Output indices;(3)Target value;(4)Weight in loss equation
dlTargetParameters = {dlTargetParameters1, dlTargetParameters2};

% Define training options
dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
dlTrainOptions('dlLambda') = 1e-2;% Fitting/Learning rate
dlTrainOptions('dlEpochs') = 1000; % Total iterations 
dlTrainOptions('dlEnhancedMomentum') = 0.25;

dlTrainOptions('dlExcludeDiverge') = 0; % Exclude non-optimals from model log
dlTrainOptions('dlCheckpointLengthCap') = 100;
dlTrainOptions('dlUpdateVerbose') = 1;
dlTrainOptions('dlTrainIncludeList') = "_IappContext";

dlTrainOptions('dlTrainRestrictList') = ["ES_noise", "INfast_noise"]; % Restrict list
dlTrainOptions('dlTrainRestrictCoef') = {.001, .001}; % Restrict coeffs
% dlTrainOptions('dlTrainSyncList') = [["ES_INfast_iGABAa_tauD", "INfast_INfast_iGABAa_tauD"]; ["ES_noise", "INfast_noise"]]; % Each array row will be overrided by synchronizing to the first element
dlTrainOptions('dlCheckpointCoefficient') = 4.0;
dlTrainOptions('dlUpdateMode') = 'batch';

dl.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

%%

dl.dlLoadOptimal 
dl.dlSimulate
dl.dlPlotAllPotentials('lfp');
dl.dlPlotBatchErrors();

%%

p = load([study_dir, '/Optimalparams.mat']);
p = p.p;
fprintf("Starting parameters > ES_noise = %f, INfast_noise = %f \n ->ES_INfast_iGABAa_tauD = %f\n->INfast_INfast_iGABAa_tauD = %f\n", 40.0, 40.0, 5.0, 5.0);
fprintf("Optimal parameters > ES_noise = %f, INfast_noise = %f \n ->ES_INfast_iGABAa_tauD = %f\n->INfast_INfast_iGABAa_tauD = %f\n", p.ES_noise, p.INfast_noise, p.ES_INfast_iGABAa_tauD, p.INfast_INfast_iGABAa_tauD);
fprintf("Optimal outputs > ES_V_fNAT = %f, ES_V_aV = %f\n", dl.dlLastOutputs{1}, dl.dlLastOutputs{2});
fprintf("Target outputs > ES_V_fNAT = %f, ES_V_aV = %f\n", target_FRQ, target_AV);

%%