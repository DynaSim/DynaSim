% DynaSim CUBA model (from Brian paper)
LIF={
  'tau=10; E=-70; thresh=-55; reset=-75; R=9; I=1.7; noise=100'
  'dV/dt=(E-V+R*I-@isyn)/tau; V(0)=-65+2*randn(1,Npop)'
  'if(V>thresh)(V=reset)'
%   'if(any(t<tspike+tabs,1))(V=reset); tabs=10'
%   'monitor V.spikes(thresh)'
};

N=8000; % Brian: 4000, 8000, 16000, 32000
time_limits=[0 1000];
compile_flag=0;

% LIF without synapses
s=[];
s.pops(1).name='LIF';
s.pops(1).size=N;
s.pops(1).equations=LIF;
data=SimulateModel(s,'time_limits',time_limits,'solver','rk1','dt',.01,'compile_flag',compile_flag,'verbose_flag',1);
figure; subplot(3,1,1); plot(data.time,data.LIF_V(:,1:min(10,N))); %ylim([-80 -50]);

%% With synapses

CLOCKampa={
  'gSYN=.1; ESYN=0; tauD=2; tauR=0.4; delay=15'
  'netcon=ones(N_pre,N_post)'
  'ISYN(X,s)=gSYN.*(s*netcon).*(X-ESYN)'
  'ds/dt=-s./tauD+((1-s)/tauR).*(1+tanh(X_pre(t-delay)/10)); s(0)=.1'
  '@isyn += -ISYN(X_post,s)'
};
EVENTampa={
  'gSYN=.1; ESYN=0; tauD=2; tauR=0.4; delay=15'
  'netcon=ones(N_pre,N_post)'
  'f(x) = (exp(-x/tauD)-exp(-x/tauR)).*(x>0)'
  'Isyn(X) = gSYN.*(sum(f(t-tspike_pre-delay))*netcon).*(X-ESYN)'
  '@isyn += Isyn(X_post)'
};
EVENTampa2={
  'gSYN=.1; ESYN=0; tauD=2; tauR=0.4; delay=15; max_int=5*tauD'
  'netcon=ones(N_pre,N_post)'
  'f(x) = (exp(-x/tauD)-exp(-x/tauR)).*(x>0)'
  'if(any(t<tspike_pre+max_int,1))(Isyn(X)=gSYN.*(sum(f(t-tspike_pre-delay))*netcon).*(X-ESYN))(Isyn(X)=0)'
  '@isyn += Isyn(X_post)'
};

netcon=zeros(N,N); % ones(N,N), zeros(N,N)

compile_flag=0;
sparse_flag=0;

% LIF with clock-based synapses
s.cons(1).direction='LIF->LIF';
s.cons(1).mechanism_list='iampa';
s.cons(1).parameters={'gSYN',.1/N,'netcon',netcon};
s.mechanisms(1).name='iampa';
s.mechanisms(1).equations=CLOCKampa;
data=SimulateModel(s,'time_limits',time_limits,'solver','rk1','dt',.01,'compile_flag',compile_flag,'sparse_flag',sparse_flag,'verbose_flag',1);
subplot(3,1,2); plot(data.time,data.LIF_V(:,1:min(10,N))); %ylim([-80 -50]);

% LIF with event-based synapses
s.mechanisms(1).equations=EVENTampa;
data=SimulateModel(s,'time_limits',time_limits,'solver','rk1','dt',.01,'compile_flag',compile_flag,'sparse_flag',sparse_flag,'verbose_flag',1);
subplot(3,1,3); plot(data.time,data.LIF_V(:,1:min(10,N))); %ylim([-80 -50]);

% LIF with event-based synapses and event-based updating
s.mechanisms(1).equations=EVENTampa2;
data=SimulateModel(s,'time_limits',time_limits,'solver','rk1','dt',.01,'compile_flag',compile_flag,'sparse_flag',sparse_flag,'verbose_flag',1);
subplot(3,1,3); plot(data.time,data.LIF_V(:,1:min(10,N))); %ylim([-80 -50]);

% N   T   compile sparse  synapse |   time (sec)
% ----------------------------------------------
% 100 1s  0       0(all)  none        0.57
% 100 1s  0       0(all)  clock       1.72
% 100 1s  0       0(all)  event(2,dt) 2.66
% ----------------------------------------------
% 100 1s  0       1(all)  none        0.57
% 100 1s  0       1(all)  clock       5.42
% 100 1s  0       1(all)  event(2,dt) 7.64
% 100 1s  0       1(none) none        0.55
% 100 1s  0       1(none) clock       4.46
% 100 1s  0       1(none) event(2,dt) 6.62
% 100 1s  1       0(all)  none        0.37
% 100 1s  1       0(all)  clock       2.55
% 100 1s  1       0(all)  event(2,dt) 3.04




% used in Brian benchmark:
% - sparse connectivity matrix
% - compilation
