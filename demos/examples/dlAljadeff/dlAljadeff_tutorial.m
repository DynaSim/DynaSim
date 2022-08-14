%% DynaSim model with local plasticity rules from Aljadeff et al., arXiv 2019
%% 
% This tutorial demonstrates how to implement the local plasticity rules introduced 
% in <https://arxiv.org/abs/1911.00307 Aljadeff et al., arXiv 2019>.
%% 
% * *Goals*: 
%% 
% # Illustrate the creation and usage of DynaSim models endowed with local plasticity 
% updates via DynaLearn. 
% # Apply the local (Hebbian, neuromodulatory and inhibitory) plasticity rules 
% introduced in Aljadeff et al., arXiv 2019 in a conductance-based Hodgkin-Huxley 
% cortical model.
%% 
% * *Structure*: After a brief introduction, the first part of the tutorial 
% shows how to load DynaSim and set up the DynaSim model. The execution of the 
% model is done in section Run and visualize the model. The second part of the 
% tutorial shows how to set up DynaLearn to dynamically adapt the weights according 
% to the local plasticity rules, run the simulation and visualize the changes. 
% Finally, two follow-up ideas are suggested for the motivated audience.
% * *Reference:*
%% 
% # Johnatan Aljadeff, James D'amour, Rachel E. Field, Robert C. Froemke, Claudia 
% Clopath. _Cortical credit assignment by Hebbian, neuromodulatory and inhibitory 
% plasticity_. arXiv:1911.00307, 2019
%% Background
% We have implemented the local plasticity rules introduced in Aljadeff et al., 
% arXiv 2019, however our aim is not to replicate the model in Aljadeff et al, 
% arXiv 2019 (see Follow-up ideas). Rather, we have preferred to illustrate how 
% these (or other) local plasticity rules can be easily implemented and used in 
% all sort of DynaSym models. Consequently, we have considered that a conductance-based 
% Hodgkin-Huxley cortical model is a great candidate to illustrate this point.
%% Load DynaSim
% If you haven't loaded DynaSim, let's do it now.

DynaSim_path = '../../..';
addpath(genpath(DynaSim_path));
close all;
clear;
clc;
%% 1. Implementation of the DynaSim model
% In this section, we specify the initial settings both of the simulation (how 
% and for how long we integrate the model) and that of the model itself (all its 
% mechanisms and their parameters). Then we will execute the model and visualize 
% its dynamics before using DynaLearn to modify its weights.
% DynaSim Parameters
% Here, we specify some parameters of the model and simulation.

% parameters
ton         = 0;
toff        = 2000;
tspan       = [ton toff]; % ms

Ne          = 40;
Ni          = 10;

g_poisson   = 2.8e-3;

tauR_GABA   = 1; % ms, rise time constant of inhibition
tauD_GABA   = 10; % ms, decay time constant of inhibition

tauR_AMPA   = 0.2; % ms, rise time constant of excitation
tauD_AMPA   = 2; % ms, decay time constant of excitation

g_rec       = 0.2;
gAMPA_ee    = 2*g_rec/Ne; % E->E
gAMPA_ei    = 1.4*g_rec/Ne; % E->I
gGABAa_ie   = 2.8*g_rec/Ni; % I->E
gGABAa_ii   = 1*g_rec/Ni; % I->I

Y_ACh       = 1; % 0; % global pairing (supposed to be trial dependent)
rho_ACh     = 0.5;
A_ACh       = 0.227;

rho_NE      = 0.012;
A_NE        = 0.5;
% Model specification
% Now, we detail the specification of the DynaSim model. Note the input mechanisms 
% for ACh and NE in E-cells. We initiliaze the weights of plastic synapses (E->E 
% and I->E) to 0.5, so they can rise and decay.

% specification

spec = [];

% dynamics
eqns = 'dV/dt = @current/Cm; Cm=1; V(0) = -90 + 20*randn(1, Npop); monitor @current, @iE, @iI';

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

spec.connections(1).direction = 'E->E'; % Later on modulated by ACh, NE, Hebbian plasticity
spec.connections(1).mechanism_list = {'iAMPA'};
spec.connections(1).parameters = {'gAMPA',gAMPA_ee,'tauR',tauR_AMPA,'tauD',tauD_AMPA,'netcon',0.5*ones(Ne,Ne)};

spec.connections(2).direction = 'E->I'; % fix
spec.connections(2).mechanism_list = {'iAMPA'};
spec.connections(2).parameters = {'gAMPA',gAMPA_ei,'tauR',tauR_AMPA,'tauD',tauD_AMPA,'netcon',ones(Ne,Ni)};

spec.connections(3).direction = 'I->E'; % Later on modulated by Inhibitory plasticity
spec.connections(3).mechanism_list = {'iGABAa'};
spec.connections(3).parameters = {'gGABAa',gGABAa_ie,'tauR',tauR_GABA,'tauD',tauD_GABA,'netcon',0.5*ones(Ni,Ne)};

spec.connections(4).direction = 'I->I'; % fix
spec.connections(4).mechanism_list = {'iGABAa'};
spec.connections(4).parameters = {'gGABAa',gGABAa_ii,'tauR',tauR_GABA,'tauD',tauD_GABA,'netcon',ones(Ni,Ni)};
% Run and visualize the model
% Finally in this section, let's run the model and visualize the dynamics of 
% the model before letting the weights change.

% Control simulation
data = dsSimulate(spec,'time_limits',tspan,'dt',0.01,'solver','euler','downsample_factor',10,'verbose_flag',1);

raster{1} = computeRaster(data.time,data.E_V);
raster{2} = computeRaster(data.time,data.I_V);

tl          = tspan;
uscaling    = 1e3; % unit scaling from kernel regression from kHz to Hz
kwidth      = 100; % width of kernel regression in ms
pool        = [Ne Ni];
time        = data.time;
dt          = time(2)-time(1);
Ts          = 1; % subsampling period in ms
flag_interp = 0;
kernel      = 'L'; % 'E'; % 'G'; % for Laplacian, Epanechnikov (default) and Gaussian
ifr         = [];
if ~all(cellfun(@isempty,raster))
    raster = raster(~cellfun(@isempty,raster));
    rate = plotRaster(tl,raster);
    for ipop = 1:numel(raster)
        [ifr(:,ipop), t] = NWKraster(time, raster{ipop}, 1:pool(ipop), kwidth, Ts, flag_interp, kernel);
    end
    ifr = uscaling*ifr;

    lineWidth   = 1;
    fontSize    = 16;
    colors      = [33, 113, 181; 239, 59, 44]./255;

    figure
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
%% 2. Create a DynaLearn object
% First, we create and initialize a DynaLearn object. Note that you can use 
% MATLAB's Coder shall your simulations are long or you plan to run many epochs 
% (see the nEpochs parameter below), but note that generating MEX files may take 
% long time.

coderFlag = false; % true; 
existFlag = false; % true;

if existFlag
    m = m.dlLoad('.');
else
    if coderFlag
        m = DynaLearn(spec, '.'); % using Matlab Coder
    else
        m = DynaLearn(spec, '.', 'raw'); % not using Matlab Coder
    end
    m.dlPathToFile = '.';
    Params = containers.Map();
    Params('tspan') = tspan;
    m.dlUpdateParams(Params);
    m.dlSave();
    system('mv dlTempFunc.m solve');
    system('mv dlTempFuncParamsChanger.m solve');
    system('mv dlFile.mat solve');
    system('mv params.mat solve');
end
% DynaLearn Parameters
% Here, we specify parameters for DynaLearn's local plasticity rules: ACh, NE, 
% Hebbian (applied to E->E) and Inhibitory (applied to I->E).

% parameters of local plasticity rules

% AMPA
iAMPA_w_min                         = 0;
iAMPA_w_max                         = 1;
fr_norm                             = 100; % used to renormalize activity [0 100] sp/s in [0 1]
f_ref                               = 0.2;

% ACh
y_ref_ACh                           = 0.005;
alpha_ACh                           = 0.575;
beta_ACh                            = 0.33;

% NE
y_ref_NE                            = 0.005;
alpha_NE                            = 0.772;

% Hebbian
alpha_Hebbian                       = 0.16;

% GABAa
iGABAa_w_min                        = 0;
iGABAa_w_max                        = 1;

% Inhibitory
a_Inhibitory                        = 0.5;
b_Inhibitory                        = 0.05;
alpha_Inhibitory                    = 0.1;

% LocalParams structure for dynalearn

Local_LR_EE_Params                  = [];
Local_LR_EE_Params.learningRules    = {'ACh','NE','Hebbian'};
Local_LR_EE_Params.source           = 'E';
Local_LR_EE_Params.target           = 'E';
Local_LR_EE_Params.connection_type  = 'iAMPA';
Local_LR_EE_Params.fr_norm          = fr_norm;
Local_LR_EE_Params.f_ref            = f_ref;
Local_LR_EE_Params.w_min            = iAMPA_w_min;
Local_LR_EE_Params.w_max            = iAMPA_w_max;
Local_LR_EE_Params.voltage          = 'V';

% ACh
Local_LR_EE_Params.y_ref_ACh        = y_ref_ACh;
Local_LR_EE_Params.Y_ACh            = Y_ACh;
Local_LR_EE_Params.rho_ACh          = rho_ACh;
Local_LR_EE_Params.alpha_ACh        = alpha_ACh;
Local_LR_EE_Params.beta_ACh         = beta_ACh;
Local_LR_EE_Params.eta_ACh          = 'iACh_eta_ACh';

% NE
Local_LR_EE_Params.y_ref_NE         = y_ref_NE;
Local_LR_EE_Params.rho_NE           = rho_NE;
Local_LR_EE_Params.alpha_NE         = alpha_NE;
Local_LR_EE_Params.eta_NE           = 'iNE_eta_NE';

% Hebbian
Local_LR_EE_Params.alpha_Hebbian    = alpha_Hebbian;
Local_LR_EE_Params.y_ref_Hebbian    = f_ref;

% Inhibitory
Local_LR_IE_Params                  = [];
Local_LR_IE_Params.learningRules    = {'Inhibitory'};
Local_LR_IE_Params.source           = 'I';
Local_LR_IE_Params.target           = 'E';
Local_LR_IE_Params.connection_type  = 'iGABAa';
Local_LR_IE_Params.fr_norm          = fr_norm;
Local_LR_IE_Params.w_min            = iGABAa_w_min;
Local_LR_IE_Params.w_max            = iGABAa_w_max;
Local_LR_IE_Params.voltage          = 'V';
Local_LR_IE_Params.iE               = 'iE';
Local_LR_IE_Params.iI               = 'iI';

Local_LR_IE_Params.a_Inhibitory     = a_Inhibitory;
Local_LR_IE_Params.b_Inhibitory     = b_Inhibitory;
Local_LR_IE_Params.alpha_Inhibitory = alpha_Inhibitory;
% TrainModel specification
% Now, we detail the specification of the DynaSim model. Note the input mechanisms 
% for ACh and NE in E-cells. We initiliaze the weights of plastic synapses (E->E 
% and I->E) to 0.5, so they can rise and decay.

dlTrainOptions                      = containers.Map();

dlTrainOptions('Local_LR_Params')   = {Local_LR_EE_Params, Local_LR_IE_Params};

nEpochs                             = 10;
dlTrainOptions('dlEpochs')          = nEpochs;
dlTrainOptions('dlBatchs')          = 1; % >1 when an error signal is defined and there are different types of trials

dlTrainOptions('dlCheckpoint')      = 'false'; % useful if an error signal is defined (it does not apply to local plasticity rules)
dlTrainOptions('dlUpdateMode')      = 'batch';
dlTrainOptions('dlLearningRule')    = 'none'; % not using non-local learning rules here
dlTrainOptions('dlSimulationFlag')  = 1; % Executing simulations (1 or 0; 1 is default, 0 only for dev debugging)
dlTrainOptions('dlOutputLogFlag')   = 1; % Autosaving trial outputs (1 or 0; 0 is default, for the case many simulations are executed)

% Not used here
dlOutputParameters                  = [];
dlInputParameters                   = {};
dlTargetParameters                  = {};

% DynaLearn train
m.dlTrain(dlInputParameters, dlOutputParameters, dlTargetParameters, dlTrainOptions);
% Run and visualize the model (revisited)
% Finally in this section, we visualize the dynamics of the model after the 
% synaptic changes. Note that we have updated the weights 10 times (nEpochs = 
% 10). Increase the number of epochs if you want to infer how the dynamics evolves 
% futher than that.

% if unsure about the index, check m.dlWeightsVariables
weights_evolution_EE = m.dlWeightsValues{1};
weights_evolution_IE = m.dlWeightsValues{2};

spec.connections(1).parameters = {'gAMPA',gAMPA_ee,'tauR',tauR_AMPA,'tauD',tauD_AMPA,'netcon',weights_evolution_EE(:,:,end)};
spec.connections(3).parameters = {'gGABAa',gGABAa_ie,'tauR',tauR_GABA,'tauD',tauD_GABA,'netcon',weights_evolution_IE(:,:,end)};

% After training simulation
data = dsSimulate(spec,'time_limits',tspan,'dt',0.01,'solver','euler','downsample_factor',10,'verbose_flag',1);

raster{1} = computeRaster(data.time,data.E_V);
raster{2} = computeRaster(data.time,data.I_V);

ifr = [];
if ~all(cellfun(@isempty,raster))
    raster = raster(~cellfun(@isempty,raster));
    rate = plotRaster(tl,raster);
    for ipop = 1:numel(raster)
        [ifr(:,ipop), t] = NWKraster(time, raster{ipop}, 1:pool(ipop), kwidth, Ts, flag_interp, kernel);
    end
    ifr = uscaling*ifr;

    figure
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
% Synaptic change
% Let's now visualize the (mean) synaptic evolution introduced by the local 
% plasticity rules.

avg_weight_evol_EE = nan(1+nEpochs,1);
avg_weight_evol_IE = nan(1+nEpochs,1);

for epoch = 0:nEpochs
    weights_epoch_EE = weights_evolution_EE(:,:,1+epoch);
    weights_epoch_IE = weights_evolution_IE(:,:,1+epoch);
    avg_weight_evol_EE(1+epoch) = mean(weights_epoch_EE(:));
    avg_weight_evol_IE(1+epoch) = mean(weights_epoch_IE(:));
end

figure
hold on
set(gca,'layer','top')
plot([0 nEpochs], [0.5 0.5],'--','color',[0.3 0.3 0.3],'linewidth',lineWidth)
plot(0:nEpochs,avg_weight_evol_EE,'color',colors(1,:),'linewidth',lineWidth)
plot(0:nEpochs,avg_weight_evol_IE,'color',colors(2,:),'linewidth',lineWidth)
xlim([0 nEpochs])
ylim([0 1])
set(gca,'fontSize',fontSize,'LineWidth',lineWidth,'TickDir','out','Box','on','XTick',0:nEpochs);
xlabel('Epochs','fontSize',fontSize)
ylabel('Avg. synaptic weight','fontSize',fontSize)
legend(' w_0',' E→E',' I→E','box','off')
%% Follow-up ideas
% We have seen how DynaLearn allows embedding local plasticity rules into arbitrary 
% DynaSim models. In fact, we have implemented the plasticity rules of <https://arxiv.org/abs/1911.00307 
% Aljadeff et al., arXiv 2019> in a conductance-based Hodgkin-Huxley cortical 
% model. From here, we have visualized two immediate follow-up projects that we 
% propose for the motivated audience: 
%% 
% # Implement the model in Aljadeff et al., arXiv 2019 in DynaSim+DynaLearn.
% # Study whether ACh, NE, Hebbian and Inhibitory plasticity rules with the 
% proper coordination are also able to approximate the Delta rule in conductance-based 
% Hodgkin-Huxley cortical models