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

% Get ready...

demos_path = findDemosPath;

% Set path to your copy of the DynaSim toolbox
dynasim_path = fullefile(demos_path, '..');

% add DynaSim toolbox to Matlab path
addpath(genpath(dynasim_path)); % comment this out if already in path

% Set where to save outputs
output_directory = fullfile(demos_path, 'outputs');
% move to root directory where outputs will be saved
cd(output_directory);

% Here we go!

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINING AND SIMULATING MODELS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Lorenz equations with phase plot

% DynaSim makes it easy to simulate arbitrary systems of ordinary
% differential equations. Simply write out the system in a cell array of
% strings, separating equations into different strings or the same string
% separated by semicolons.

eqns={
  's=10; r=27; b=2.666';
  'dx/dt=s*(y-x)';
  'dy/dt=r*x-y-x*z';
  'dz/dt=-b*z+x*y';
};
data=SimulateModel(eqns,'tspan',[0 100],'ic',[1 2 .5],'solver','rk4');
% tspan: time limits on integration [ms]
% ic: initial conditions
% solver: numerical method to use (default: rk4 = "4th-order Runge-Kutta")

% All models are numerically integrated using a DynaSim solver function
% created uniquely for a given model and stored in a directory named
% "solve". The file that solves the system (i.e,. numerically integrates
% it) is stored in data.simulator_options and can be viewed or rerun afterwards:
edit(data.simulator_options.solve_file)

% Every component of the model is assigned to a "population", and the 
% population name (default: 'pop1') is prepended to all variable and
% function names.

% Simulated data can be easily plotted using the resulting data structure:
figure; plot(data.pop1_x,data.pop1_z); 
title('Lorenz equations'); xlabel('x'); ylabel('z')

%% Izhikevich neuron with noisy drive 
% (reference: p274 of "Dynamical Systems in Neuroscience" by Izhikevich)

% The DynaSim data structure always contains the model state variables,
% time vector, and a copy of the DynaSim model structure that was
% simulated. Additionally, functions can be recorded and returned in the
% DynaSim data structure if indicated using the "monitor" keyword.
% Syntax: monitor FUNCTION

eqns={
  'C=100; vr=-60; vt=-40; k=.7; Iapp=70; ton=200; toff=800';
  'a=.03; b=-2; c=-50; d=100; vpeak=35';
  'dv/dt=(k*(v-vr)*(v-vt)-u+I(t))/C; v(0)=vr';
  'du/dt=a*(b*(v-vr)-u); u(0)=0';
  'if(v>vpeak)(v=c; u=u+d)';
  'I(t)=Iapp*(t>ton&t<toff)*(1+.5*rand)'; % define applied input using reserved variables 't' for time and 'dt' for fixed time step of numerical integration
  'monitor I';                            % indicate to store applied input during simulation
};
data=SimulateModel(eqns,'tspan',[0 1000]);
% plot the simulated voltage and monitored input function
figure; 
subplot(2,1,1); plot(data.time,data.pop1_v); % plot voltage
xlabel('time (ms)'); ylabel('v'); title('Izhikevich neuron')
subplot(2,1,2); plot(data.time,data.pop1_I); % plot input function
xlabel('time (ms)'); ylabel('Iapp');

% note: "t", "dt", and "T" are special variables that can be used in model
% equations. "t" represents the current time point of the simulation. 
% "dt" is the fixed time step used for numeric integration. "T" is the full
% simulated time vector defined before simulation begins.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUNNING SETS OF SIMULATIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 'vary' indicates the variable to vary, the values it should take, and the 
% object (population or connection) whose variable should be varied. 

% Syntax 1: vary={{object, variable, value1},{object, variable, value2},...}
%   - this is useful for simulating an arbitrary set of parameter values
% Syntax 2: vary={object, variable, values; ...}
%   - this is useful for varying parameters systematically (described later)

% Izhikevich study of neuro-computational properties (using Syntax 1)
% based on: http://www.izhikevich.org/publications/izhikevich.m
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
data=SimulateModel(eqns,'tspan',[0 250],'vary',vary);
PlotData(data);

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
data=SimulateModel(eqns);
figure; plot(data.time,data.(data.labels{1}))
xlabel('time (ms)'); ylabel('membrane potential (mV)'); title('Hodgkin-Huxley neuron')

% Equivalent Hodgkin-Huxley neuron with predefined mechanisms
data=SimulateModel('dv/dt=10+@current/Cm; Cm=1; v(0)=-65; {iNa,iK}');
figure; plot(data.time,data.(data.labels{1}))
xlabel('time (ms)'); ylabel('membrane potential (mV)'); title('Hodgkin-Huxley neuron')

% View the mechanism files:
[~,eqnfile]=LocateModelFiles('iNa.mech'); edit(eqnfile{1});
[~,eqnfile]=LocateModelFiles('iK.mech');  edit(eqnfile{1});
% Mechanisms can be custom built; however, DynaSim does come pakaged with
% some common ones like popular ion currents (see <dynasim>/models).

% Example of a bursting neuron model using three active current mechanisms:
eqns='dv/dt=5+@current; {iNaF,iKDR,iM}; gNaF=100; gKDR=5; gM=1.5; v(0)=-70';
data=SimulateModel(eqns,'tspan',[0 200]);
figure; plot(data.time,data.(data.labels{1}))
xlabel('time (ms)'); ylabel('membrane potential (mV)'); title('Intrinsically Bursting neuron')

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
data=SimulateModel(s);
PlotData(data);
PlotData(data,'variable',{'E_v','E_I_iGABAa_ISYN'});

% View the connection mechanism file:
[~,eqnfile]=LocateModelFiles('iAMPA.mech'); edit(eqnfile{1});

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVING SIMULATED DATA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How to: set 'save_data_flag' to 1 and (optionally) 'study_dir' to /path/to/outputs

%% Save data from a single simulation
% Example using the previous sPING model:
data=SimulateModel(s,'save_data_flag',1,'study_dir','demo_sPING_1');

%% Save data from a set of simulations

% Specify what to vary 
% Tip: use 'vary' Syntax 2 to systematically vary a parameter
vary={'E','Iapp',[0 10 20]}; % vary the amplitude of tonic input to E-cells
data=SimulateModel(s,'save_data_flag',1,'study_dir','demo_sPING_2',...
                     'vary',vary);

% load and plot the saved data
data_from_disk = ImportData('demo_sPING_2');
PlotData(data_from_disk);
PlotData(data_from_disk,'variable','E_v');

% Vary a connection parameter
vary={'I->E','tauD',[5 10 15]}; % inhibition decay time constant from I to E

% Vary two parameters (run a simulation for all combinations of values)
vary={
  'E'   ,'Iapp',[0 10 20];      % amplitude of tonic input to E-cells
  'I->E','tauD',[5 10 15]       % inhibition decay time constant from I to E
  };
SimulateModel(s,'save_data_flag',1,'study_dir','demo_sPING_3',...
                'vary',vary,'verbose_flag',1);
data=ImportData('demo_sPING_3');
PlotData(data);
PlotData(data,'plot_type','rastergram');
PlotData(data,'plot_type','power');
PlotFR(data); % examine how mean firing rate changes with Iapp and tauD

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUNNING SIMULATIONS ON THE CLUSTER
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How to: set 'cluster_flag' to 1
% Requirement: you must be logged on to a cluster that recognizes 'qsub'

% Run three simulations in parallel jobs and save the simulated data
eqns='dv/dt=@current+I; {iNa,iK}';
vary={'','I',[0 10 20]};
SimulateModel(eqns,'save_data_flag',1,'study_dir','demo_cluster_1',...
                   'vary',vary,'cluster_flag',1,'overwrite_flag',1,'verbose_flag',1);
% tips for checking job status:
% !qstat -u <YOUR_USERNAME>
% !cat ~/batchdirs/demo_cluster_1/pbsout/sim_job1.out
data=ImportData('demo_cluster_1');
PlotData(data);

% Repeat but also save plotted data
eqns='dv/dt=@current+I; {iNa,iK}';
vary={'','I',[0 10 20]};
SimulateModel(eqns,'save_data_flag',1,'study_dir','demo_cluster_2',...
                   'vary',vary,'cluster_flag',1,'overwrite_flag',1,'verbose_flag',1,...
                   'plot_functions',@PlotData);
% !cat ~/batchdirs/demo_cluster_2/pbsout/sim_job1.out

% Save multiple plots and pass custom options to each plotting function
eqns='dv/dt=@current+I; {iNa,iK}';
vary={'','I',[0 10 20]};
SimulateModel(eqns,'save_data_flag',1,'study_dir','demo_cluster_3',...
                   'vary',vary,'cluster_flag',1,'overwrite_flag',1,'verbose_flag',1,...
                   'plot_functions',{@PlotData,@PlotData},...
                   'plot_options',{{},{'plot_type','power'}});
% !cat ~/batchdirs/demo_cluster_3/pbsout/sim_job1.out
 
% Post-simulation analyses can be performed similarly by passing
% analysis function handles and options using 'analysis_functions' and
% 'analysis_options'. 

% Note: options will be passed to plot and analysis functions in the order
% given. You can pass handles and options for any built-in, pre-packaged,
% or custom functions.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MORE FEATURES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulating large models can be sped up significantly by compiling the
% simulation before running it. DynaSim makes this easy to do using the
% 'compile_flag' option in SimulateModel. Note: compiling the model can 
% take several seconds to minutes; however, it only compiles the first time
% it is run and is significantly faster on subsequent runs.

data=SimulateModel(s,'compile_flag',1);
PlotData(data);

% Now run again:
data=SimulateModel(s,'compile_flag',1);
PlotData(data);

