% DynaLearn applied to PFC cell
% addpath(genpath('/Users/jason/Documents/me/docs - research/andre/BetaGammaPushPull/DynaSim'));

% Optimize RMP by varying Eleak

RunDuration = 500; % In ms (miliseconds)
study_dir = 'dlModels/HH';

spec = dsCheckSpecification('HH');
vary = {'HH','Eleak',[-70 -60 -50 -40 -30]};
solver_options = {'tspan',[0 RunDuration],'solver','rk1','dt',.01,'compile_flag',0,'verbose_flag',1};
data = dsSimulate(spec,'vary',vary,solver_options{:});
dsPlot(data,'plot_type','waveform');

target_RMP = -62;

dl = DynaLearn(spec, study_dir, 'raw', 'HH');
dl.dlSave();

% Define optimization parameters
dlInputParameters = {dlNullInputs(RunDuration)}; % Null inputs (no external stimulation).
dlOutputParameters = [{['HH', '_V'], 1, [100 RunDuration], 'av'}]; % 'av' is average voltage
dlTargetParameters = {[{'MQE', 1, target_RMP, 1}]}; % Format: (1)Mode;(2)Output indices;(3)Target value;(4)Weight in loss equation

% Define training options
dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
dlTrainOptions('dlLambda') = 1e-3; % Fitting/Learning rate
dlTrainOptions('dlEpochs') = 100; % Total iterations 
dlTrainOptions('dlEnhancedMomentum') = 0.3;
dlTrainOptions('dlTrainExcludeList') = {'_Cm', '_Iapp', '_gNa', '_gK', '_gleak', '_ENa', '_EK'};
dlTrainOptions('dlCheckpointCoefficient') = 1.1;
dl.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

dl.dlPlotBatchErrors(1)
dl.dlLastErrorsLog
dl.dlLoadOptimal
dl.dlSimulate
dl.dlPlotAllPotentials('lfp')

optimal_param_file = fullfile(study_dir, 'Optimalparams.mat');
load(optimal_param_file, 'p'); p
fprintf('Optimal Eleak: %g\n', p.HH_ileak_Eleak)

% Initial parameters:
%               HH_Cm: 1
%             HH_Iapp: 0
%          HH_iNa_gNa: 120
%          HH_iNa_ENa: 50
%     HH_iNa_IC_noise: 0
%         HH_iNa_h_IC: 0
%         HH_iNa_m_IC: 0
%            HH_iK_gK: 36
%            HH_iK_EK: -77
%      HH_iK_IC_noise: 0
%          HH_iK_n_IC: 0
%      HH_ileak_gleak: 0.3000
%      HH_ileak_Eleak: -70
%             HH_Npop: 1

% Optimal parameters:
%                 HH_Cm: 1
%               HH_Iapp: 0
%            HH_iNa_gNa: 120
%            HH_iNa_ENa: 50
%       HH_iNa_IC_noise: 0
%           HH_iNa_h_IC: 0
%           HH_iNa_m_IC: 0
%              HH_iK_gK: 36
%              HH_iK_EK: -77
%        HH_iK_IC_noise: 0
%            HH_iK_n_IC: 0
%        HH_ileak_gleak: 0.3000
%        HH_ileak_Eleak: -46.0587
%               HH_Npop: 1

%% Optimizing average firing rate (FR) by varying noise level

RunDuration = 500; % In ms (miliseconds)
study_dir = 'dlModels/noisyHH';

% add noise mechanisms
spec = dsApplyModifications('HH.pop',{'HH','mechanism_list','+noise'});
vary = {'HH','amp',[-100 100 200 300 400 500 1000]}; % % 20Hz is around amp=250
solver_options = {'tspan',[0 RunDuration],'solver','rk1','dt',.01,'compile_flag',0,'verbose_flag',1};
data = dsSimulate(spec,'vary',vary,solver_options{:});
dsPlotFR(data); % ~ 80 spk/s
% data.model.parameters

%% set baseline noise level
spec = dsApplyModifications(spec,{'HH','noise_amp', 1000}); 
 
% optimize FR
target_FR = 20; % Hz
dl = DynaLearn(spec, study_dir, 'raw', 'HH');
dl.dlSave();

% Define optimization parameters
dlInputParameters = {dlNullInputs(RunDuration)}; % Null inputs (no external stimulation).
dlOutputParameters = [{['HH', '_V'], 1, [100 RunDuration], 'afr'}; {['HH', '_V'], 1, [100 RunDuration], 'av'}];
dlTargetParameters = {[{'MSE', 1, target_FR, 1}; {'MSE', 2, -62, 1}]}; % Format: (1)Mode;(2)Output indices;(3)Target value;(4)Weight in loss equation

% Define training options
dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
dlTrainOptions('dlLambda') = 1e-2; % Fitting/Learning rate
dlTrainOptions('dlEpochs') = 100; % Total iterations 
dlTrainOptions('dlEnhancedMomentum') = 0.5;

dlTrainOptions('dlExcludeDiverge') = 0; % Exclude non-optimals from model log
dlTrainOptions('dlCheckpointLengthCap') = 11;
% TODO: !dlTrainIncludeList
dlTrainOptions('dlTrainExcludeList') = {'_Cm', '_Iapp', '_gNa', '_gK', '_gleak', ...
  '_ENa', '_EK', '_Eleak', '_IC_noise'};
dlTrainOptions('dlCheckpointCoefficient') = 1.7;
dl.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

dl.dlLoadOptimal 
dl.dlSimulate
try dl.dlPlotAllPotentials('lfp'); end

optimal_spec = load(fullfile(study_dir, 'Optimalobject.mat'));
optimal_param_file = fullfile(study_dir, 'Optimalparams.mat');
load(optimal_param_file, 'p'); p
fprintf('Optimal noise: %g\n', p.HH_noise_noise_amp)

% dl.dlResetTraining

% todo: 
 % - dl.dlLoadOptimal: load optimal dsData, dlModel structures & params, replace optimal object
 % - replace o1 = obj.dlApplyIFRKernel(dlOutput); with o1 = obj.dlMeanFR(dlOutput);


% return

%% Optimizing average firing rate (FR) by varying Iapp

RunDuration = 500; % In ms (miliseconds)
study_dir = 'dlModels/iHH';

% add noise mechanisms
spec = dsCheckSpecification('HH');
vary = {'HH','Iapp',[-10 -4 -1 0 1 4 7 8 8.80720 8.80721 8.80741 9 10 14 17 20]};
solver_options = {'tspan',[0 RunDuration],'solver','rk1','dt',.01,'compile_flag',0,'verbose_flag',1};
data = dsSimulate(spec,'vary',vary,solver_options{:});
dsPlotFR(data);
% data.model.parameters

% spec = dsApplyModifications(spec,{'HH','noise_amp', 1000}); 
 
% optimize FR
target_FR = 80; % Hz
dl = DynaLearn(spec, study_dir, 'raw', 'HH');
dl.dlSave();

%%
% Define optimization parameters
dlInputParameters = {dlNullInputs(RunDuration)}; % Null inputs (no external stimulation).
dlOutputParameters = [{['HH', '_V'], 1, [100 RunDuration], 'afr'}; {['HH', '_V'], 1, [100 RunDuration], 'av'}];
dlTargetParameters = {[{'MSE', 1, target_FR, 1}; {'MSE', 2, -55, 1}]}; % Format: (1)Mode;(2)Output indices;(3)Target value;(4)Weight in loss equation

% Define training options
dlTrainOptions = containers.Map(); % Train options; MUST be a map data structure
dlTrainOptions('dlLambda') = 1e-3; % Fitting/Learning rate
dlTrainOptions('dlEpochs') = 100; % Total iterations 
dlTrainOptions('dlEnhancedMomentum') = 0.25;

dlTrainOptions('dlExcludeDiverge') = 0; % Exclude non-optimals from model log
dlTrainOptions('dlCheckpointLengthCap') = 5;
dlTrainOptions('dlUpdateVerbose') = 1;
dlTrainOptions('dlTrainIncludeList') = ["Iapp", "_Cm", "noise"];

dlTrainOptions('dlTrainExcludeList') = ["_Cm", "noise", "_gNa", "_gK", "_gleak", ...
  "_ENa", "_EK", "_Eleak", "_IC_noise"];
dlTrainOptions('dlCheckpointCoefficient') = 1.74;
dl.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);

%%
dl.dlLoadOptimal 
dl.dlSimulate
try dl.dlPlotAllPotentials('lfp'); end

optimal_spec = load(fullfile(study_dir, 'Optimalobject.mat'));
optimal_param_file = fullfile(study_dir, 'Optimalparams.mat');
load(optimal_param_file, 'p'); p
fprintf('Optimal noise: %g\n', p.HH_noise_noise_amp)

% dl.dlResetTraining

% todo: 
 % - dl.dlLoadOptimal: load optimal dsData, dlModel structures & params, replace optimal object
 % - replace o1 = obj.dlApplyIFRKernel(dlOutput); with o1 = obj.dlMeanFR(dlOutput);


return



%%

Ne=1;       % number of E-cells per layer
Ni=0;  % number of I-cells per layer (FS + RSNP)
spec = get_PFC_1layer('DS02PYjs',Ne,'DS02FSjs',Ni/2,'DS02RSNPjs',0);
solver_options={'tspan',[0 500],'solver','rk1','dt',.01,'compile_flag',0,'verbose_flag',1};
vary = [];
data=dsSimulate(spec,'vary',vary,solver_options{:});
dsPlot(data,'plot_type','waveform');



% % Add Poisson-based inputs to superficial PY dendrite (20kHz w/gAMPA=1e-3uS)
% generic input and state equations
input_def={'input(V)=iAMPA(V); monitor input; onset=50; offset=inf;';
           'iAMPA(V)=-gAMPA.*sAMPA(k,:).*(V-EAMPA); EAMPA=0; gAMPA=0;';
           'sAMPA=getPoissonGating(0,dcAMPA,acAMPA,freq,0,onset,offset,tauAMPA,T,Npop); dcAMPA=0; acAMPA=0; tauAMPA=2; freq=0';
          }; % {'gAMPA',1e-3 to 1e-5,'dcAMPA',20e3,'Iapp',.1}; % [gAMPA]=uS, [dcAMPA]=Hz
state_equations=['dV/dt=(@current+input(V)+Iapp*(t>onset&t<offset))./Cm; Cm=1; Iapp=0; V(0)=-65;' input_def{:}];
spec = dsApplyModifications(spec,{'Ed','equations',state_equations});

solver_options={'tspan',[0 1500],'solver','rk1','dt',.01,'compile_flag',0,'verbose_flag',1};
tic
vary = {'Ed','dcAMPA',[0 10]*1e3; 'Ed','gH',[0, 50, 500];'Ed','gAMPA',1e-4};
data=dsSimulate(spec,'vary',vary,solver_options{:});
toc % 
dsPlot(data,'plot_type','waveform');
tic; dsPlot(data,'plot_type','power'); toc
