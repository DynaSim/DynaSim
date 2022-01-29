% Single layer PFC network with PY, FS, and RSNP cells
% Based on a NEURON model converted to a DynaSim implementation

% Single layer PFC model size
Ne=8;  % # pyramidal cells
Nfs=2; % # fast-spiking (FS) cells
Nrs=2; % # regular-spiking non-pyramidal (RSNP) cells

% -------------------------------------------------------------------------
% Specify populations: Custom cell models (see also: 'help get_PFC_cell')
% -------------------------------------------------------------------------
% Two-compartment Pyramidal cell model ('Es'=soma, 'Ed'=dendrite)
spec=get_PFC_cell('DS02PYjs',Ne);
% One-compartment FS model ('FS' = PV+ interneuron inhibiting PY soma)
spec.populations(end+1)=getfield(get_PFC_cell('DS02FSjs',Nfs),'populations');
% One-compartment RSNP model ('RSNP' = CB+ interneuron inhibiting PY dendrite)
spec.populations(end+1)=getfield(get_PFC_cell('DS02RSNPjs',Nrs),'populations');

% -------------------------------------------------------------------------
% Specify network connections (all-to-all)
% -------------------------------------------------------------------------
% [DS02,DG07]: AMPA (taur=.2,taud=1,E=0), NMDA (taur=2.3,taud=95,E=0), GABA (taur=.5,taud=5,E=-75)
tauAMPAr=.2;  % ms, AMPA rise time
tauAMPAd=1;   % ms, AMPA decay time
tauNMDAr=2.3; % ms, NMDA rise time
tauNMDAd=95;  % ms, NMDA decay time
tauGABAr=.5;  % ms, GABAa rise time
tauGABAd=5;   % ms, GABAa decay time
% [DS02]:
Npc=100;  % # principal cells in original [DS02] publication
Nin=37;   % # interneurons in original [DS02] publication
gAMPAee=3e-3;       % uS, PY->PY, maximal AMPA conductance
gNMDAee=gAMPAee/50; % uS, PY->PY, maximal NMDA conductance
gGABAie=.2e-3;      % uS, FS->PY, maximal GABAa conductance
gAMPAei=.74e-3;     % uS, PY->IN
gNMDAei=gAMPAei/50; % uS, PY->IN
gGABAii=.6e-3;      % uS, FS->IN

% [D97] DeFelipe, J. (1997). Types of neurons, synaptic connections and chemical characteristics of cells immunoreactive for calbindin-D28K, parvalbumin and calretinin in the neocortex. Journal of chemical neuroanatomy, 14(1), 1-19.
% - NMDA colocalizes with PV+ but not CB+
% - PV+ INs target PY soma; CB+ INs target PY dendrite
% - PV+ INs are FS; CB+ INs are RSNP or LTS
% [Traub/Whittington] (2010, in "Cortical Oscillations in Health and Disease")
% - CB+ inhibition is longer-lasting than PV+ inhibition
tauGABAdCB=13; % ms, GABAa decay time for inhibition from RSNP CB+ interneurons

% recurrent connections between pyramidal cells (add to existing intercompartmental connections)
index=find(strcmp('Es->Ed',{spec.connections.direction}),1,'first');
spec.connections(index).mechanism_list={'iAMPA','iNMDA',spec.connections(index).mechanism_list{:}};
spec.connections(index).parameters={'gAMPA',gAMPAee*Npc/Ne,'gNMDA',gNMDAee*Npc/Ne,'EAMPA',0,'ENMDA',0,...
  'tauAMPAr',tauAMPAr,'tauAMPA',tauAMPAd,'tauNMDAr',tauNMDAr,'tauNMDA',tauNMDAd,spec.connections(index).parameters{:}};
% pyramidal<->FS connections
spec.connections(end+1).direction='Es->FS';
spec.connections(end).mechanism_list={'iAMPA','iNMDA'};
spec.connections(end).parameters={'gAMPA',gAMPAei*Npc/Ne,'gNMDA',gNMDAei*Npc/Ne,'EAMPA',0,'ENMDA',0,...
  'tauAMPAr',tauAMPAr,'tauAMPA',tauAMPAd,'tauNMDAr',tauNMDAr,'tauNMDA',tauNMDAd};
spec.connections(end+1).direction='FS->Es';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABAie*Nin/Nfs,'tauGABAr',tauGABAr,'tauGABA',tauGABAd,'EGABA',-75};
% pyramidal<->RSNP connections
spec.connections(end+1).direction='Es->RSNP';
spec.connections(end).mechanism_list={'iAMPA'};
spec.connections(end).parameters={'gAMPA',gAMPAei*Npc/Ne,'EAMPA',0,'tauAMPAr',tauAMPAr,'tauAMPA',tauAMPAd};
spec.connections(end+1).direction='RSNP->Ed';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABAie*Nin/Nrs,'tauGABAr',tauGABAr,'tauGABA',tauGABAdCB,'EGABA',-75};
% interneuron<->interneuron connections
spec.connections(end+1).direction='FS->FS';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABAii*Nin/Nfs,'tauGABAr',tauGABAr,'tauGABA',tauGABAd,'EGABA',-75};
spec.connections(end+1).direction='RSNP->RSNP';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABAii*Nin/Nrs,'tauGABAr',tauGABAr,'tauGABA',tauGABAdCB,'EGABA',-75};
spec.connections(end+1).direction='FS->RSNP';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABAii*Nin/Nfs,'tauGABAr',tauGABAr,'tauGABA',tauGABAd,'EGABA',-75};
spec.connections(end+1).direction='RSNP->FS';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABAii*Nin/Nrs,'tauGABAr',tauGABAr,'tauGABA',tauGABAdCB,'EGABA',-75};

%% update inputs (background poisson + rhythmic population signal)

% randomize seed for Poisson input
rng('shuffle'); 
seed=getfield(rng,'Seed');

TargetPop='Ed'; % name of population(s) to stimulate (for inputs and recurrent connections Es->Ed)
OutputPop='Es'; % name of population to analyze (for post-simulation analysis)
EiPop='Es';     % excitatory population input compartment for E/I connectivity
EoPop='Es';     % excitatory population output compartment for E/I connectivity
IPop='FS';      % inhibitory population input=output compartment for E/I connectivity

TargetPopIndex=find(strcmp(TargetPop,{spec.populations.name}));
OutputPopIndex=find(strcmp(OutputPop,{spec.populations.name}));
IPopIndex=find(strcmp(IPop,{spec.populations.name}));

% generic inputs and dependent state equations
input_def={'input(V)=iNoise(V)+iSignal(V)+iPoisson(V); monitor input,ginput,iNoise,iSignal,iPoisson; onset=400; offset=inf;';
           'ginput(t)=gNoise.*sNoise(k,:)+gSignal.*sSignal(k,:)+gPoisson.*sPoisson(k,:);';
           'iNoise(V)=gNoise.*sNoise(k,:).*(V-0); gNoise=.001; dcNoise=5000;';
           'iSignal(V)=gSignal.*sSignal(k,:).*(V-0); gSignal=.02;';
           'iPoisson(V)=gPoisson.*sPoisson(k,:).*(V-0); gPoisson=0;';           
           'sNoise=getPopRhythmGating(1,Npop,0,dcNoise,0,0,T,0,inf,2,1,seed);';
           'sSignal=getPopRhythmGating(Ninputs,Npop,pSignal,rSignal,fSignal,wSignal,T,onset,offset,tauSignal,1,seed,ones(1,Npop),actype);';
           'sPoisson=getPoissonGating(0,rSignal*Ninputs,rSignal*Ninputs,fSignal,0,onset,offset,tauSignal,T,Npop,ones(1,Npop));';
           'rSignal=5; wSignal=5; fSignal=32; pSignal=1; tauSignal=2; Ninputs=100; seed=0; actype=''step'';';
           };
state_equations=['dV/dt=(@current-input(V)+Iapp*(t>onset&t<offset))./Cm; Cm=1; Iapp=0; V(0)=-65;' input_def{:}];

% common signal parameters
Ninputs=100; % # of simulated input cells
onset=400;
offset=inf;
% E-cell inputs
% noise:
gNoiseE=.0007;
dcNoiseE=1000; % spk/s, independent homogeneous Poisson (e.g., 1000 cells at 5 spk/s)
% signal:
gSignalE=.0015;
pSignalE=0;  % probability that input i connects to target j (set to 0 for independent inputs)
wSignalE=1; % ms, synchrony of population spikes on a given cycle
fSignalE=25; % Hz, population frequency
rSignalE=10; % spk/s, average cell firing rate
tauSignalE=2; % ms, synaptic decay constant
% FS-cell inputs
% noise:
gNoiseI=.001;
dcNoiseI=0; % spk/s, independent homogeneous Poisson (e.g., 1000 cells at 5 spk/s)
% signal:
gSignalI=0;
pSignalI=0;  % probability that input i connects to target j (set to 0 for independent inputs)
wSignalI=1; % ms, synchrony of population spikes on a given cycle
fSignalI=25; % Hz, population frequency
rSignalI=0; % spk/s, average cell firing rate
tauSignalI=2; % ms, synaptic decay constant

% collect input and biophysical parameters
E_input_parameters ={'Ninputs',Ninputs,'onset',onset,'offset',offset,'gSignal',gSignalE,'rSignal',rSignalE,'wSignal',wSignalE,'fSignal',fSignalE,'pSignal',pSignalE,'tauSignal',tauSignalE,'gNoise',gNoiseE,'dcNoise',dcNoiseE,'gPoisson',0,'seed',seed};
I_input_parameters ={'Ninputs',Ninputs,'onset',onset,'offset',offset,'gSignal',gSignalI,'rSignal',rSignalI,'wSignal',wSignalI,'fSignal',fSignalI,'pSignal',pSignalI,'tauSignal',tauSignalI,'gNoise',gNoiseI,'dcNoise',dcNoiseI,'gPoisson',0,'seed',seed};

% update population specification (add inputs and input parameters)
spec=dsApplyModifications(spec,{TargetPop,'equations',state_equations}); % Ed
spec.populations(TargetPopIndex).parameters=cat(2,E_input_parameters,dsRemoveKeyval(spec.populations(TargetPopIndex).parameters,E_input_parameters(1:2:end)));
spec=dsApplyModifications(spec,{IPop,'equations',state_equations}); % FS
spec.populations(IPopIndex).parameters=cat(2,I_input_parameters,dsRemoveKeyval(spec.populations(IPopIndex).parameters,I_input_parameters(1:2:end)));

% -------------------------------------------------------------------------
% Run simulations
% -------------------------------------------------------------------------

% drive network with different frequencies
vary = {'Ed','fSignal',0:5:50};
tspan = [0 1000];  % [beg end], ms
dt = .01;         % fixed time step, ms
compile_flag = 0; % whether to compile simulation
verbose_flag = 1; % whether to display verbose logs
solver_options = {'tspan',tspan,'solver','rk1','dt',dt,'compile_flag',compile_flag,'verbose_flag',verbose_flag};
data = dsSimulate(spec,'vary',vary,solver_options{:});

% plot waveforms
dsPlot(data,'plot_type','waveform');

% plot mean firing rate vs. varied parameter (input frequency)
select_data = dsSelect(data,'time_limits',[onset, tspan(2)]);
dsPlotFR(select_data);
