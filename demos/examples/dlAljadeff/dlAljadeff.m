%% DynaLearn model using learning rules from Aljadeff et al., arXiv 2019
%% Reference: Aljadeff et al. Cortical credit assignment by Hebbian, neuromodulatory and inhibitory plasticity. arXiv:1911.00307, 2019

close all;
clear;
clc;

% dynasim

% path to dynasim
dynasim_path = '/home/jsardid/Documents/tr3ball/code/mcode/DynaSim';
addpath(genpath(dynasim_path));

% parameters
ton = 0;
toff = 2000;
tspan = [ton toff]; % ms
transient = 1000;

Ne = 40; % 8;
Ni = 10; % 2;

g_poisson = 2.8e-3; % 3e-3;

tauR_GABA = 1; % ms, rise time constant of inhibition
tauD_GABA = 10; % ms, decay time constant of inhibition

tauR_AMPA = 0.2; % ms, rise time constant of excitation
tauD_AMPA = 2; % ms, decay time constant of excitation

g_rec = 0.4;
gAMPA_ee = 1*g_rec/Ne; % 1.7*g_rec/Ne; % E->E within layer
gAMPA_ei = 1.4*g_rec/Ne; % 1.4*g_rec/Ne; % E->I within layer
gGABAa_ie = 1.4*g_rec/Ni; % 2*g_rec/Ni; % I->E within layer
gGABAa_ii = 1*g_rec/Ni; % 1.4*g_rec/Ni; % I->I within layer

Y_ACh = 1; % 0; % global pairing
rho_ACh = 0.5; % 0.825;
A_ACh = 0.227;

rho_NE = 0.012;
A_NE = 0.5; % 1.605;

% dynamics
eqns = 'dV/dt = @current/Cm; Cm=1; V(0) = -90 + 20*randn(1, Npop); monitor @current, @iE, @iI';

% specification

spec = [];

% populations

% E-cells
spec.populations(1).name = 'E';
spec.populations(1).size = Ne;
spec.populations(1).equations = eqns;
spec.populations(1).mechanism_list = {'iL','iNa','iK','iPoisson','iACh','iNE'};
spec.populations(1).parameters = {'gext',g_poisson,'Y_ACh',Y_ACh,'rho_ACh',rho_ACh,'A_ACh',A_ACh,'rho_NE',rho_NE,'A_NE',A_NE};

% I-cells
spec.populations(2).name = 'I';
spec.populations(2).size = Ni;
spec.populations(2).equations = eqns;
spec.populations(2).mechanism_list = {'iL','iNa','iK','iPoisson'};
spec.populations(2).parameters = {'gext',g_poisson};

% connections

spec.connections(1).direction = 'E->E'; % ACh, NE, Hebbian plasticity
spec.connections(1).mechanism_list = {'iAMPA'};
spec.connections(1).parameters = {'gAMPA',gAMPA_ee,'tauR',tauR_AMPA,'tauD',tauD_AMPA,'netcon',0.5*ones(Ne,Ne)};

spec.connections(2).direction = 'E->I'; % fix
spec.connections(2).mechanism_list = {'iAMPA'};
spec.connections(2).parameters = {'gAMPA',gAMPA_ei,'tauR',tauR_AMPA,'tauD',tauD_AMPA,'netcon',0.5*ones(Ne,Ni)};

spec.connections(3).direction = 'I->E'; % Inhibitory plasticity
spec.connections(3).mechanism_list = {'iGABAa'};
spec.connections(3).parameters = {'gGABAa',gGABAa_ie,'tauR',tauR_GABA,'tauD',tauD_GABA,'netcon',0.5*ones(Ni,Ne)};

spec.connections(4).direction = 'I->I'; % fix
spec.connections(4).mechanism_list = {'iGABAa'};
spec.connections(4).parameters = {'gGABAa',gGABAa_ii,'tauR',tauR_GABA,'tauD',tauD_GABA,'netcon',0.5*ones(Ni,Ni)};

% Control simulation
data = dsSimulate(spec,'time_limits',tspan,'dt',0.01,'solver','euler','downsample_factor',10,'verbose_flag',1); % ,'mex_flag',1);

raster{1} = computeRaster(data.time,data.E_V);
raster{2} = computeRaster(data.time,data.I_V);

% raster{1} = computeRasterSorted(data.time,data.E_V);
% raster{2} = computeRasterSorted(data.time,data.I_V);

tl = tspan; % [transient, tspan(2)];
uscaling = 1e3; % unit scaling from kernel regression from kHz to Hz
kwidth = 500; % width of kernel regression in ms
pool = [Ne Ni];
time = data.time;
dt = time(2)-time(1);
Ts = 1; % double(dt); % 1;  % subsampling period in ms
flag_interp = 0;
kernel = 'L'; % 'E'; % 'G';
time = data.time;
ifr = [];
if ~all(cellfun(@isempty,raster))
    raster = raster(~cellfun(@isempty,raster));
    rate = plotRaster(tl,raster);
    for ipop = 1:numel(raster)
        [ifr(:,ipop), t] = NWKraster(time, raster{ipop}, 1:pool(ipop), kwidth, Ts, flag_interp, kernel);
    end
    ifr = uscaling*ifr;

    lineWidth = 1;
    fontSize = 16;
    colors = [33, 113, 181; 239, 59, 44]./255;

    figure('visible','on')
    hold on
    set(gca,'layer','top')
    for ipop = 1:numel(raster)
        plot(t,ifr(:,ipop),'color',colors(ipop,:),'linewidth',lineWidth)
    end
    xlim([0 max(t)])
    set(gca,'fontSize',fontSize,'LineWidth',lineWidth,'TickDir','out','Box','on','YTick',0:max([4,200/kwidth]):max(ifr(:)));
    xlabel('Time (ms)','fontSize',fontSize)
    ylabel('Inst. firing rate (sp/s)','fontSize',fontSize)
else
    disp('empty rasters')
end

%%% dynalearn

coderFlag = false; % true;
existFlag = false; % true;

if existFlag
    %% Loading the DynaLearn Class (if it already exists)
    m = m.dlLoad('.'); % ~ 10sec
    m.dlStudyDir = 'solve';
else
    %% Creating the DynaLearn Class
    if coderFlag
        m = DynaLearn(spec, '.'); % using Matlab Coder (~ hour)
    else
        m = DynaLearn(spec, '.', 'raw'); % not using Matlab Coder (~ minute)
    end
    m.dlPathToFile = '.';
    Params = containers.Map();
    Params('tspan') = tspan;
    m.dlUpdateParams(Params);
    m.dlSave(); % (< second)
    system('mv dlTempFunc.m solve');
    system('mv dlTempFuncParamsChanger.m solve');
    system('mv dlFile.mat solve');
    system('mv params.mat solve');
end

%% Simulation and general plotting
% m.dlSimulate(); % (optional) simulate it (~ minute)
% m.dlPlotAllPotentials('ifr'); % Instantaneous firing rate
% m.dlPlotAllPotentials('lfp'); % Local field potential

% parameters of local plasticity rules

% AMPA
iAMPA_w_min = 0;
iAMPA_w_max = 1;

fr_norm = 100; % used to renormalize [0 100] sp/s in [0 1]
f_ref = 0.2;

% ACh
y_ref_ACh = 0.01;
alpha_ACh = 0.575;
beta_ACh = 0.331;

% NE
y_ref_NE = 0.01;
alpha_NE = 0.772;

% Hebbian
alpha_Hebbian = 0.016;

% GABAa
iGABAa_w_min = 0;
iGABAa_w_max = 1;

% Inhibitory
a_Inhibitory = 1;
b_Inhibitory = 0.05;
alpha_Inhibitory = 0.638;

% LocalParams structure for dynalearn

Local_LR_EE_Params = [];
Local_LR_EE_Params.learningRules = {'ACh','NE','Hebbian'};
Local_LR_EE_Params.source = 'E';
Local_LR_EE_Params.target = 'E';
Local_LR_EE_Params.connection_type = 'iAMPA';
Local_LR_EE_Params.fr_norm = fr_norm;
Local_LR_EE_Params.f_ref = f_ref;
Local_LR_EE_Params.w_min = iAMPA_w_min;
Local_LR_EE_Params.w_max = iAMPA_w_max;
Local_LR_EE_Params.voltage = 'V';
Local_LR_EE_Params.spikes = 'V';

% ACh
Local_LR_EE_Params.y_ref_ACh = y_ref_ACh;
Local_LR_EE_Params.Y_ACh = Y_ACh;
Local_LR_EE_Params.rho_ACh = rho_ACh;
Local_LR_EE_Params.alpha_ACh = alpha_ACh;
Local_LR_EE_Params.beta_ACh = beta_ACh;
Local_LR_EE_Params.eta_ACh = 'iACh_eta_ACh';

% NE
Local_LR_EE_Params.y_ref_NE = y_ref_NE;
Local_LR_EE_Params.rho_NE = rho_NE;
Local_LR_EE_Params.alpha_NE = alpha_NE;
Local_LR_EE_Params.eta_NE = 'iNE_eta_NE';

% Hebbian
Local_LR_EE_Params.alpha_Hebbian = alpha_Hebbian;
Local_LR_EE_Params.y_ref_Hebbian = f_ref;

% Inhibitory
Local_LR_IE_Params = [];
Local_LR_IE_Params.learningRules = {'Inhibitory'};
Local_LR_IE_Params.source = 'I';
Local_LR_IE_Params.target = 'E';
Local_LR_IE_Params.connection_type = 'iGABAa';
Local_LR_IE_Params.fr_norm = fr_norm;
Local_LR_IE_Params.w_min = iGABAa_w_min;
Local_LR_IE_Params.w_max = iGABAa_w_max;
Local_LR_IE_Params.voltage = 'V';
Local_LR_IE_Params.iE = 'iE';
Local_LR_IE_Params.iI = 'iI';

Local_LR_IE_Params.a_Inhibitory = a_Inhibitory;
Local_LR_IE_Params.b_Inhibitory = b_Inhibitory;
Local_LR_IE_Params.alpha_Inhibitory = alpha_Inhibitory;

dlTrainOptions = containers.Map();

dlTrainOptions('Local_LR_Params') = {Local_LR_EE_Params, Local_LR_IE_Params};

dlTrainOptions('dlEpochs') = 2; % 100;
dlTrainOptions('dlBatchs') = 1;

dlTrainOptions('dlCheckpoint') = 'false';
dlTrainOptions('dlUpdateMode') = 'batch';
dlTrainOptions('dlLearningRule') = 'none'; % not using non-local learning rules % 'BioDeltaRule'; % DeltaRule, BioDeltaRule, RWDelta, ...

dlTrainOptions('dlSimulationFlag') = 1; % Manually turning simulation, on or off (on is default and recommended)
dlTrainOptions('dlOutputLogFlag') = 1; % Autosaving trial outputs, on or off (off is default and recommended) % TODO Output/Random/SameValueProblem

%% Train
dlOutputParameters = [];
dlInputParameters = {};
dlTargetParameters = {};

m.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);
% m.dlPlotAllPotentials('ifr');

%% Run a simulation (after training)

% m.dlResetTraining(); % Reset logs and optimal state error (but not the optimal state file)
% m.dlLoadOptimal();  % Load the current optimal state (if exists)
% m.dlRunSimulation();
% m.dlPlotAllPotentials('ifr');
