%% Make DynaSim Autogen Files

warning('off','MATLAB:lang:cannotClearExecutingFunction');

cwd = pwd;

% Set where to save outputs
output_directory = fullfile(dsGetConfig('ds_data_path'), 'autoGenSave_temp');

% Erase output_directory
if isdir(output_directory)
  rmdir(output_directory, 's')
end

% move to root directory where outputs will be saved
mkdirSilent(output_directory);
cd(output_directory);

autogenOptions = {'auto_gen_test_data_flag',1, 'random_seed', 1, 'verbose_flag',0, 'visible', 'off'};

%% Lorenz equations
eqns={
  's=10; r=27; b=2.666';
  'dx/dt=s*(y-x)';
  'dy/dt=r*x-y-x*z';
  'dz/dt=-b*z+x*y';
};
data=dsSimulate(eqns, 'tspan',[0 100], 'ic',[1 2 .5], 'solver','rk4', 'study_dir','demo_lorenz', autogenOptions{:});
dsPlot(data, autogenOptions{:}); close gcf

%% Izhikevich neuron with noisy drive
eqns={
  'C=100; vr=-60; vt=-40; k=.7; Iapp=70; ton=10; toff=20';
  'a=.03; b=-2; c=-50; d=100; vpeak=35';
  'dv/dt=(k*(v-vr)*(v-vt)-u+I(t))/C; v(0)=vr';
  'du/dt=a*(b*(v-vr)-u); u(0)=0';
  'if(v>vpeak)(v=c; u=u+d)';
  'I(t)=Iapp*(t>ton&t<toff)*(1+.5*rand)'; % define applied input using reserved variables 't' for time and 'dt' for fixed time step of numerical integration
  'monitor I';                            % indicate to store applied input during simulation
};
data=dsSimulate(eqns, 'tspan',[0 100], 'study_dir','demo_izhikevich', autogenOptions{:});
dsPlot(data, autogenOptions{:}); close gcf

%% RUNNING SETS OF Izhikevich SIMULATIONS
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

%% normal
data=dsSimulate(eqns, 'tspan',[0 100], 'vary',vary, 'study_dir','demo_izhikevich_vary', autogenOptions{:});

% TODO: fix dsPlot bug
% dsPlot(data, autogenOptions{:}); close gcf

%% parfor
data=dsSimulate(eqns, 'tspan',[0 100], 'vary',vary, 'study_dir','demo_izhikevich_vary_parfor',...
  'parallel_flag',1, autogenOptions{:});

%% parfor compile
data=dsSimulate(eqns, 'tspan',[0 100], 'vary',vary, 'study_dir','demo_izhikevich_vary_parfor_compile',...
  'parallel_flag',1, 'compile_flag',1, autogenOptions{:});

%% compile
data=dsSimulate(eqns, 'tspan',[0 100], 'vary',vary, 'study_dir','demo_izhikevich_vary_compile',...
  'compile_flag',1, autogenOptions{:});

%% Hodgkin-Huxley neuron equations (without predefined mechanisms)
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
data=dsSimulate(eqns, 'study_dir','demo_hh_1', autogenOptions{:});

dsPlot(data, autogenOptions{:}); close gcf

% Equivalent Hodgkin-Huxley neuron with predefined mechanisms
data=dsSimulate('dv/dt=10+@current/Cm; Cm=1; v(0)=-65; {iNa,iK}', 'study_dir','demo_hh_2', autogenOptions{:});

% Example of a bursting neuron model using three active current mechanisms:
eqns='dv/dt=5+@current; {iNaF,iKDR,iM}; gNaF=100; gKDR=5; gM=1.5; v(0)=-70';
data=dsSimulate(eqns, 'tspan',[0 100], 'study_dir','demo_hh_3', autogenOptions{:});
dsPlot(data, autogenOptions{:}); close gcf

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
s.populations(1).size=8;
s.populations(1).equations=eqns;
s.populations(1).mechanism_list={'iNa','iK'};
s.populations(1).parameters={'Iapp',5,'gNa',120,'gK',36,'noise',40};
s.populations(2).name='I';
s.populations(2).size=2;
s.populations(2).equations=eqns;
s.populations(2).mechanism_list={'iNa','iK'};
s.populations(2).parameters={'Iapp',0,'gNa',120,'gK',36,'noise',40};
s.connections(1).direction='I->E';
s.connections(1).mechanism_list={'iGABAa'};
s.connections(1).parameters={'tauD',10,'gSYN',.1,'netcon','ones(N_pre,N_post)'};
s.connections(2).direction='E->I';
s.connections(2).mechanism_list={'iAMPA'};
s.connections(2).parameters={'tauD',2,'gSYN',.1,'netcon',ones(8,2)};

data=dsSimulate(s, 'study_dir','demo_sPING_0', autogenOptions{:});
dsPlot(data, autogenOptions{:}); close all

%% Covary
% Vary two parameters (run a simulation for all combinations of values)
vary={'(E,I)','(EK1,EK2)',[-80 -60]};
%       This sets modifications:
%           * E_EK1, E_EK2, I_EK1, I_EK2 = -80
%           * E_EK1, E_EK2, I_EK1, I_EK2 = -60
data = dsSimulate(s, 'study_dir','demo_sPING_covary1', 'vary',vary, autogenOptions{:});
              
vary={'(E,I)','(EK1,EK2)',[-80 -60; -85 -65]};
%       This sets modifications:
%           * E_EK1, I_EK1 = -80 and E_EK2, I_EK2 = -85
%           * E_EK1, I_EK1 = -60 and E_EK2, I_EK2 = -65
data = dsSimulate(s, 'study_dir','demo_sPING_covary2', 'vary',vary, autogenOptions{:});
              
vary={'(E,I)','(EK1,EK2)',cat(3,[-80 -60], [-85 -65])};
%       This sets modifications:
%           * E_EK1, E_EK2 = -80 and I_EK1, I_EK2 = -85
%           * E_EK1, E_EK2 = -60 and I_EK1, I_EK2 = -65
data = dsSimulate(s,'study_dir','demo_sPING_covary3', 'vary',vary, autogenOptions{:});
              
vary={'(E,I)','(EK1,EK2)',cat(3, [-75 -55; -80 -60], [-85 -65; -90 -70])};
%       This sets modifications:
%           * E_EK1 = -75, E_EK2 = -80, I_EK1 = -85, I_EK2 = -90.
%           * E_EK1 = -55, E_EK2 = -60, I_EK1 = -65, I_EK2 = -70.
data = dsSimulate(s,'study_dir','demo_sPING_covary4', 'vary',vary, autogenOptions{:});


%% SAVING SIMULATED DATA

%% Save data from a single simulation
% Example using the previous sPING model:
data=dsSimulate(s, 'save_data_flag',1, 'study_dir','demo_sPING_1', autogenOptions{:});

%% Save data from a set of simulations

% Specify what to vary
% Tip: use 'vary' Syntax 2 to systematically vary a parameter
vary={'E','Iapp',[0 10 20]}; % vary the amplitude of tonic input to E-cells
data=dsSimulate(s, 'save_data_flag',1, 'study_dir','demo_sPING_2',...
                     'vary',vary, autogenOptions{:});

% load and plot the saved data
% data_from_disk = dsImport('demo_sPING_2', 'auto_gen_test_data_flag',1);
% dsPlot(data_from_disk, autogenOptions{:}); close all
% dsPlot(data_from_disk,'variable','E_v', autogenOptions{:}); close all

% Vary two parameters (run a simulation for all combinations of values)
vary={
  'E'   ,'Iapp',[0 10 20];      % amplitude of tonic input to E-cells
  'I->E','tauD',[5 10 15]       % inhibition decay time constant from I to E
  };
data=dsSimulate(s, 'save_data_flag',1, 'study_dir','demo_sPING_3',...
                'vary',vary, 'verbose_flag',1, autogenOptions{:});
% data=dsImport('demo_sPING_3', autogenOptions{:});
dsPlot(data, autogenOptions{:}); close all
dsPlot(data,'plot_type','rastergram', autogenOptions{:}); close all
dsPlot(data,'plot_type','power', autogenOptions{:}); close all
dsPlotFR(data, autogenOptions{:}); close all

%% RUNNING SIMULATIONS ON THE CLUSTER

% Run three simulations in parallel jobs and save the simulated data
eqns='dv/dt=@current+I; {iNa,iK}';
vary={'','I',[0 10 20]};
dsSimulate(eqns,'save_data_flag',1,'study_dir','demo_cluster_1_loop',...
                   'vary',vary,'cluster_flag',1, autogenOptions{:});

% Array mode
eqns='dv/dt=@current+I; {iNa,iK}';
vary={'','I',[0 10 20]};
dsSimulate(eqns,'save_data_flag',1,'study_dir','demo_cluster_1_array',...
                   'vary',vary,'cluster_flag',1, 'qsub_mode', 'array',...
                   autogenOptions{:});

% Repeat but also save plotted data
eqns='dv/dt=@current+I; {iNa,iK}';
vary={'','I',[0 10 20]};
dsSimulate(eqns,'save_data_flag',1,'study_dir','demo_cluster_2',...
                   'vary',vary,'cluster_flag',1, autogenOptions{:},...
                   'plot_functions',@dsPlot);

% Save multiple plots and pass custom options to each plotting function
eqns='dv/dt=@current+I; {iNa,iK}';
vary={'','I',[0 10 20]};
dsSimulate(eqns,'save_data_flag',1,'study_dir','demo_cluster_3',...
                   'vary',vary,'cluster_flag',1, autogenOptions{:},...
                   'plot_functions',{@dsPlot,@dsPlot},...
                   'plot_options',{{},{'plot_type','power'}});

%% Compilation

% Compile
data=dsSimulate(s,'compile_flag',1, 'study_dir','demo_sPING_3_compile', autogenOptions{:});

% Combine compilation and parallelization
vary={'','I',0:2:14};
data=dsSimulate(s, 'compile_flag',1, 'parallel_flag',1, 'vary', vary, 'study_dir','demo_sPING_3_compile_parallel', autogenOptions{:});
dsPlot(data, autogenOptions{:}); close all

%% Matlab Solvers
% Run three simulations in parallel jobs and save the simulated data
eqns='dv/dt=@current+I; {iNa,iK}';
vary={'','I',[0 10 20]};

% od45
data = dsSimulate(eqns,'save_data_flag',1,'study_dir','demo_cluster_1_ml_solver_ode45',...
                   'vary',vary,'cluster_flag',0,'overwrite_flag',1,...
                   'dt',0.01, 'downsample_factor',100,'solver','ode45',...
                    'matlab_solver_options', {'InitialStep', 0.01}, 'compile_flag',0, autogenOptions{:});
  
% ode45 compiled
data = dsSimulate(eqns,'save_data_flag',1,'study_dir','demo_cluster_1_ml_solver_ode45_compiled',...
                   'vary',vary,'cluster_flag',0,'overwrite_flag',1,...
                   'dt',0.01, 'downsample_factor',100,'solver','ode45',...
                    'matlab_solver_options', {'InitialStep', 0.01}, 'compile_flag',1, autogenOptions{:});
               
% ode23
data = dsSimulate(eqns,'save_data_flag',1,'study_dir','demo_cluster_1_ml_solver_ode23s',...
                   'vary',vary,'cluster_flag',0,'overwrite_flag',1,...
                   'dt',0.01, 'downsample_factor',100,'solver','ode23s',...
                    'matlab_solver_options', {'InitialStep', 0.01}, 'compile_flag',0, autogenOptions{:});
                  
% ode23 compiled
data = dsSimulate(eqns,'save_data_flag',1,'study_dir','demo_cluster_1_ml_solver_ode23s_compiled',...
                   'vary',vary,'cluster_flag',0,'overwrite_flag',1,...
                   'dt',0.01, 'downsample_factor',100,'solver','ode23s',...
                    'matlab_solver_options', {'InitialStep', 0.01}, 'compile_flag',1, autogenOptions{:});

%% one file mode
study_dir=fullfile('.', 'study_HH_varyI_one_file');
eqns='dv/dt=@current+I; {iNa,iK}';
vary={'','I',[0:10:30]};

dsSimulate(eqns,'vary',vary, 'study_dir',study_dir,'save_data_flag',1,...
  'cluster_flag',1, 'qsub_mode','array', 'one_solve_file_flag',1,...
  'tspan',[0 100], 'dt',0.01, 'solver','euler', 'downsample_factor', 100,...
  'plot_functions',{@dsPlot}, 'plot_options',{{'format', 'jpg', 'visible', 'off'}},...
  autogenOptions{:});

%% Remove output_directory
cd(cwd)

fprintf('Removing temporary output_directory: %s\n', output_directory)
rmdir(output_directory, 's')

fprintf('\nDone making autogen data.\n\n')