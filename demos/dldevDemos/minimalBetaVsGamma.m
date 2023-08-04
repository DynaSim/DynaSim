% addpath(genpath('/Users/jason/Documents/me/docs - research/andre/BetaGammaPushPull/DynaSim'));

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=[];
s.populations(1).name='ES';
s.populations(1).size=80;
s.populations(1).equations=eqns;
s.populations(1).mechanism_list={'iNa','iK'};
s.populations(1).parameters={'Iapp',Iapp,'gNa',120,'gK',36,'noise',40};
s.populations(2).name='INfast';
s.populations(2).size=20;
s.populations(2).equations=eqns;
s.populations(2).mechanism_list={'iNa','iK'};
s.populations(2).parameters={'Iapp',0,'gNa',120,'gK',36,'noise',40};
s.connections(1).direction='INfast->ES';
s.connections(1).mechanism_list={'iGABAa'};
s.connections(1).parameters={'tauD',tauGABAfast,'gGABAa',gGABAfast};
s.connections(2).direction='ES->INfast';
s.connections(2).mechanism_list={'iAMPA'};
s.connections(2).parameters={'tauD',tauAMPA,'gAMPA',gAMPA};
s.connections(3).direction='INfast->INfast';
s.connections(3).mechanism_list={'iGABAa'};
s.connections(3).parameters={'tauD',tauGABAfast,'gGABAa',gGABAfast};

vary = {'(INfast->INfast, INfast->ES)','tauD',[1 5 10 20 50]};
tic
data = dsSimulate(s, 'tspan', [0 1000], 'downsample_factor', 10, ...
                     'solver', 'euler', 'dt', .002, ...
                     'compile_flag', 0, 'parallel_flag', 1, ...
                     'vary', vary, 'verbose_flag', 1);
toc

dsPlot(data);
dsPlot(data, 'plot_type','raster');
tic; dsPlot(data, 'plot_type','power'); toc


%% 2. Optimize contextual input for beta/gamma push-pull

% Optimize wContext to ES, INfast, INslow 
% Condition 1 (context=0): gamma>beta
% Condition 2 (context=1): beta>gamma
%   where gamma = 55-75Hz, beta = 20-40Hz

condition1('_Context') = 0; % unpredictable
condition2('_Context') = 1; % predictable
dlInputs = {condition1, condition2};
  
dlTrainOptions('dlTrainIncludeList') = {'_wContext'}; % scalar [-20,20], to ES, INfast, INslow

% targets (Rpenalty)
% condition1: gamma > beta
% condition2: beta > gamma

% define equations of cell model (same for E and I populations)
eqns={
  'dV/dt=IappCue+wContext.*Context+@current+noise*randn(1,N_pop); Context=1; IappCue=0; wContext=0; noise=0'
};

IappCue = 10;
wContext = 10;

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
s.populations(1).parameters={'IappCue',IappCue,'gNa',120,'gK',36,'noise',40};
s.populations(2).name='INfast';
s.populations(2).size=20;
s.populations(2).equations=eqns;
s.populations(2).mechanism_list={'iNa','iK'};
s.populations(2).parameters={'IappCue',IappCue,'gNa',120,'gK',36,'noise',40};
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

vary = {'INfast','IappCue',[0 IappCue]; 'INslow','wContext',[0 wContext/2 wContext]};
tic
data = dsSimulate(s, 'tspan', [0 1000], 'downsample_factor', 10, ...
                     'solver', 'euler', 'dt', .002, ...
                     'compile_flag', 0, 'parallel_flag', 1, ...
                     'vary', vary, 'verbose_flag', 1);
toc

dsPlot(data);
dsPlot(data, 'plot_type','raster');
tic; dsPlot(data, 'plot_type','power'); toc

