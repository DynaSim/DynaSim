% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DynaSim Demos
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Download DynaSim toolbox from https://github.com/dynasim/dynasim.
  or, download using git: git clone https://github.com/dynasim/dynasim.git
For further documentation, see tutorial.m in the demos directory
Sign up for user mailing list at: https://groups.google.com/forum/#!forum/dynasim-users.

Tip: In Matlab, you can obtain more information associated with any function "FUNCTION_NAME"
by entering "help FUNCTION_NAME" in the command window. Use the "See also" list
at the end of the help section to browse through related help documentation.
%}

% Add DynaSim to path if it's not already there
if exist('setupDynaSimPath','file')
    setupDynaSimPath;
else
    error('Add the DynaSim folder to the MATLAB path - e.g. run addpath(genpath(DynaSimPath))');
end

% Set where to save outputs
output_directory = ds.getConfig('demos_path');

% move to root directory where outputs will be saved
mkdirSilent(output_directory);
cd(output_directory);

% Here we go!

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINING AND SIMULATING MODELS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Lorenz equations with phase plot
eqns={
  's=10; r=27; b=2.666';
  'dx/dt=s*(y-x)';
  'dy/dt=r*x-y-x*z';
  'dz/dt=-b*z+x*y';
};
data=dsSimulate(eqns, 'tspan',[0 100], 'ic',[1 2 .5], 'solver','rk4', 'study_dir','demo_lorenz', 'auto_gen_test_data_flag',1);

%% Izhikevich neuron with noisy drive
eqns={
  'C=100; vr=-60; vt=-40; k=.7; Iapp=70; ton=200; toff=800';
  'a=.03; b=-2; c=-50; d=100; vpeak=35';
  'dv/dt=(k*(v-vr)*(v-vt)-u+I(t))/C; v(0)=vr';
  'du/dt=a*(b*(v-vr)-u); u(0)=0';
  'if(v>vpeak)(v=c; u=u+d)';
  'I(t)=Iapp*(t>ton&t<toff)*(1+.5*rand)'; % define applied input using reserved variables 't' for time and 'dt' for fixed time step of numerical integration
  'monitor I';                            % indicate to store applied input during simulation
};
data=dsSimulate(eqns, 'tspan',[0 1000], 'study_dir','demo_izhikevich', 'auto_gen_test_data_flag',1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUNNING SETS OF SIMULATIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eqns={
  'a=.02; b=.2; c=-65; d=6; I=14';
  'dv/dt=.04*v^2+5*v+140-u+I; v(0)=-70';
  'du/dt=a*(b*v-u); u(0)=-20';
  'if(v>=30)(v=c;u=u+d)';
  };
P='pop1'; % name of population
vary={
  {P,'a',.02; P,'b',.2 ; P,'c',-50; P,'d',2;  P,'I',15} % tonic bursting
  {P,'a',.01; P,'b',.2 ; P,'c',-65; P,'d',8;  P,'I',30} % spike frequency adaptation
  {P,'a',.02; P,'b',.2 ; P,'c',-65; P,'d',6;  P,'I',7}  % spike latency
  {P,'a',.03; P,'b',.25; P,'c',-52; P,'d',0;  P,'I',0}  % rebound burst
  {P,'a',1;   P,'b',1.5; P,'c',-60; P,'d',0;  P,'I',-65}% bistability
  {P,'a',.02; P,'b',1  ; P,'c',-55; P,'d',4;  P,'I',1}  % accomodation
  {P,'a',-.02;P,'b',-1 ; P,'c',-60; P,'d',8;  P,'I',80} % inhibition-induced spiking
  {P,'a',-.026;P,'b',-1; P,'c',-45; P,'d',0;  P,'I',70} % inhibition-induced bursting
  };
data=dsSimulate(eqns, 'tspan',[0 250], 'vary',vary, 'study_dir','demo_izhikevich_vary', 'auto_gen_test_data_flag',1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QUICKLY BUILDING LARGE MODELS FROM EXISTING "MECHANISMS"
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mechanisms are predefined reusable sub-models meant to be incorporated 
% in larger complete models. Examples of mechanisms in neuron models include
% ion currents, pumps, etc. Once defined, they can be easily incorporated
% into larger models by simply listing the name of the file containing their
% equations, without the need to re-write any of the mechanism equations.
% This greatly simplifies large model prototyping and re-configuration.

% Hodgkin-Huxley neuron equations (without predefined mechanisms)
eqns={
  'gNa=120; gK=36; Cm=1';
  'INa(v,m,h) = gNa.*m.^3.*h.*(v-50)';
  'IK(v,n) = gK.*n.^4.*(v+77)';
  'dv/dt = (10-INa(v,m,h)-IK(v,n))/Cm; v(0)=-65';
  'dm/dt = aM(v).*(1-m)-bM(v).*m; m(0)=.1';
  'dh/dt = aH(v).*(1-h)-bH(v).*h; h(0)=.1';
  'dn/dt = aN(v).*(1-n)-bN(v).*n; n(0)=0';
  'aM(v) = (2.5-.1*(v+65))./(exp(2.5-.1*(v+65))-1)';
  'bM(v) = 4*exp(-(v+65)/18)';
  'aH(v) = .07*exp(-(v+65)/20)';
  'bH(v) = 1./(exp(3-.1*(v+65))+1)';
  'aN(v) = (.1-.01*(v+65))./(exp(1-.1*(v+65))-1)';
  'bN(v) = .125*exp(-(v+65)/80)';
};
data=dsSimulate(eqns, 'study_dir','demo_hh_1', 'auto_gen_test_data_flag',1);

% Equivalent Hodgkin-Huxley neuron with predefined mechanisms
data=dsSimulate('dv/dt=10+@current/Cm; Cm=1; v(0)=-65; {iNa,iK}', 'study_dir','demo_hh_2', 'auto_gen_test_data_flag',1);

% Example of a bursting neuron model using three active current mechanisms:
eqns='dv/dt=5+@current; {iNaF,iKDR,iM}; gNaF=100; gKDR=5; gM=1.5; v(0)=-70';
data=dsSimulate(eqns, 'tspan',[0 200], 'study_dir','demo_hh_3', 'auto_gen_test_data_flag',1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILDING LARGE MODELS WITH MULTIPLE POPULATIONS AND CONNECTIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sparse Pyramidal-Interneuron-Network-Gamma (sPING)

% define equations of cell model (same for E and I populations)
eqns={
  'dv/dt=Iapp+@current+noise*randn(1,N_pop)';
  'monitor iGABAa.functions, iAMPA.functions'
};
% Tip: monitor all functions of a mechanism using: monitor MECHANISM.functions

% create DynaSim specification structure
s=[];
s.populations(1).name='E';
s.populations(1).size=80;
s.populations(1).equations=eqns;
s.populations(1).mechanism_list={'iNa','iK'};
s.populations(1).parameters={'Iapp',5,'gNa',120,'gK',36,'noise',40};
s.populations(2).name='I';
s.populations(2).size=20;
s.populations(2).equations=eqns;
s.populations(2).mechanism_list={'iNa','iK'};
s.populations(2).parameters={'Iapp',0,'gNa',120,'gK',36,'noise',40};
s.connections(1).direction='I->E';
s.connections(1).mechanism_list={'iGABAa'};
s.connections(1).parameters={'tauD',10,'gSYN',.1,'netcon','ones(N_pre,N_post)'};
s.connections(2).direction='E->I';
s.connections(2).mechanism_list={'iAMPA'};
s.connections(2).parameters={'tauD',2,'gSYN',.1,'netcon',ones(80,20)};
data=dsSimulate(s, 'study_dir','demo_sPING_0', 'auto_gen_test_data_flag',1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVING SIMULATED DATA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How to: set 'save_data_flag' to 1 and (optionally) 'study_dir' to /path/to/outputs

%% Save data from a single simulation
% Example using the previous sPING model:
data=dsSimulate(s, 'save_data_flag',1, 'study_dir','demo_sPING_1', 'auto_gen_test_data_flag',1);

%% Save data from a set of simulations

% Specify what to vary
% Tip: use 'vary' Syntax 2 to systematically vary a parameter
vary={'E','Iapp',[0 10 20]}; % vary the amplitude of tonic input to E-cells
data=dsSimulate(s, 'save_data_flag',1, 'study_dir','demo_sPING_2',...
                     'vary',vary, 'auto_gen_test_data_flag',1);

% load and plot the saved data
data_from_disk = ds.importData('demo_sPING_2', 'auto_gen_test_data_flag',1);
% dsPlot(data_from_disk);
% dsPlot(data_from_disk,'variable','E_v');

% Vary a connection parameter
vary={'I->E','tauD',[5 10 15]}; % inhibition decay time constant from I to E

% Vary two parameters (run a simulation for all combinations of values)
vary={
  'E'   ,'Iapp',[0 10 20];      % amplitude of tonic input to E-cells
  'I->E','tauD',[5 10 15]       % inhibition decay time constant from I to E
  };
dsSimulate(s, 'save_data_flag',1, 'study_dir','demo_sPING_3',...
                'vary',vary, 'verbose_flag',1, 'auto_gen_test_data_flag',1);
data=ds.importData('demo_sPING_3', 'auto_gen_test_data_flag',1);
% dsPlot(data);
% dsPlot(data,'plot_type','rastergram');
% dsPlot(data,'plot_type','power');
% ds.plotFR(data); % examine how mean firing rate changes with Iapp and tauD

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUNNING SIMULATIONS ON THE CLUSTER
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run three simulations in parallel jobs and save the simulated data
% eqns='dv/dt=@current+I; {iNa,iK}';
% vary={'','I',[0 10 20]};
% dsSimulate(eqns,'save_data_flag',1,'study_dir','demo_cluster_1',...
%                    'vary',vary,'cluster_flag',1,'overwrite_flag',1,'verbose_flag',1);
% % tips for checking job status:
% % !qstat -u <YOUR_USERNAME>
% % !cat ~/batchdirs/demo_cluster_1/pbsout/sim_job1.out
% data=ds.importData('demo_cluster_1');
% dsPlot(data);
% 
% % Repeat but also save plotted data
% eqns='dv/dt=@current+I; {iNa,iK}';
% vary={'','I',[0 10 20]};
% dsSimulate(eqns,'save_data_flag',1,'study_dir','demo_cluster_2',...
%                    'vary',vary,'cluster_flag',1,'overwrite_flag',1,'verbose_flag',1,...
%                    'plot_functions',@dsPlot);
% % !cat ~/batchdirs/demo_cluster_2/pbsout/sim_job1.out
% 
% % Save multiple plots and pass custom options to each plotting function
% eqns='dv/dt=@current+I; {iNa,iK}';
% vary={'','I',[0 10 20]};
% dsSimulate(eqns,'save_data_flag',1,'study_dir','demo_cluster_3',...
%                    'vary',vary,'cluster_flag',1,'overwrite_flag',1,'verbose_flag',1,...
%                    'plot_functions',{@dsPlot,@dsPlot},...
%                    'plot_options',{{},{'plot_type','power'}});

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MORE FEATURES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulating large models can be sped up significantly by compiling the
% simulation before running it. DynaSim makes this easy to do using the
% 'compile_flag' option in dsSimulate. Note: compiling the model can 
% take several seconds to minutes; however, it only compiles the first time
% it is run and is significantly faster on subsequent runs.

data=dsSimulate(s,'compile_flag',1, 'study_dir','demo_sPING_3_compile_1', 'auto_gen_test_data_flag',1);
% dsPlot(data);
