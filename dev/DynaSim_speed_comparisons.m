%% LIF population vs HH population

HH={
  'gNa=120; gK=36; Cm=1';
  'INa(v,m,h) = gNa.*m.^3.*h.*(v-50)';
  'IK(v,n) = gK.*n.^4.*(v+77)';
  'dv/dt = (I-INa(v,m,h)-IK(v,n))/Cm; v(0)=-65; I=10';
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
LIF={
  'dV/dt=(E-V+R*I+noise*randn-@isyn)/tau; V(0)=-65'
  'if(V>thresh)(V=reset)'
  'tau=10; E=-70; thresh=-55; reset=-75; R=9; I=1.55; noise=100'
  'monitor V.spikes(thresh)'
};   

Ne=2;%8000;%800;%80;%4;%2;%1
Ni=2;%2000;%200;%20;%4;%2;%1

% Single HH without synapses
s=[];
s.pops(1).name='E';
s.pops(1).size=Ne;
s.pops(1).equations=HH;
s.pops(2).name='I';
s.pops(2).size=Ni;
s.pops(2).equations=HH;
s.pops(2).parameters={'I',0,'noise',0};
data=SimulateModel(s,'time_limits',[0 200],'solver','rk1','dt',.01,'verbose_flag',1);
PlotData(data); set(gcf,'position',[0.3203 0.6241 0.2641 0.1843]);
% Ne=Ni=1: .049 sec
% Ne=Ni=2: .278 sec
% Ne=Ni=4: .304 sec
% Ne=80,Ni=20: .956 sec
% Ne=800,Ni=200: 7.682 sec
% Ne=8000,Ni=2000: 85.031 sec

% LIF without synapses
s=[];
s.pops(1).name='E';
s.pops(1).size=Ne;
s.pops(1).equations=LIF;
s.pops(2).name='I';
s.pops(2).size=Ni;
s.pops(2).equations=LIF;
s.pops(2).parameters={'I',0,'noise',0};
data=SimulateModel(s,'time_limits',[0 200],'solver','rk1','dt',.01,'verbose_flag',1);
PlotData(data); set(gcf,'position',[0.3203 0.6241 0.2641 0.1843]);
% Ne=Ni=1: .133 sec
% Ne=Ni=2: .196 sec
% Ne=Ni=4: .216 sec
% Ne=80,Ni=20: .619 sec
% Ne=800,Ni=200, record 100 spikes: 4.506 sec
% Ne=8000,Ni=2000, record 100 spikes: 45.162 sec
% Ne=800,Ni=200, record 5 spikes: .664 sec
% Ne=8000,Ni=2000, record 5 spikes: 7.258 sec


%% LIF network: clock-based vs event-based stepwise synapses with full connectivity
% i.e., do matrix multiplication (s*netcon) for all elements on every time step

CLOCKampa={
  'gSYN=0.1; ESYN=0; tauD=2; tauR=0.4'
  'netcon=ones(N_pre,N_post)'
  'ISYN(X,s)=(gSYN.*(s*netcon).*(X-ESYN))'
  'ds/dt=-s./tauD+((1-s)/tauR).*(1+tanh(X_pre/10)); s(0)=.1'
  '@isyn += -ISYN(X_post,s)'
};
EVENTampa={
  'gSYN=.5; ESYN=0; tauD=2; tauR=0.4; delay=15'
  'netcon=ones(N_pre,N_post)'
  'f(x) = (exp(-x/tauD)-exp(-x/tauR)).*(x>0)'
  'Isyn(X) = gSYN.*(sum(f(t-tspike_pre-delay))*netcon).*(X-ESYN)'
  '@isyn += Isyn(X_post)'
};
EVENTampa={
  'gSYN=.5; ESYN=0; tauD=2; tauR=0.4; delay=15'
  'netcon=ones(N_pre,N_post)'
  'f(x) = max(0,exp(-x/tauD)-exp(-x/tauR))'
  'Isyn(X) = gSYN.*(sum(f(t-tspike_pre-delay))*netcon).*(X-ESYN)'
  '@isyn += Isyn(X_post)'
};

% HH with clock-based synapses
s=[];
s.pops(1).name='E';
s.pops(1).size=Ne;
s.pops(1).equations=HH;
s.pops(2).name='I';
s.pops(2).size=Ni;
s.pops(2).equations=HH;
s.pops(2).parameters={'I',0,'noise',0};
s.cons(1).direction='E->I';
s.cons(1).mechanism_list='iampa';
s.mechanisms(1).name='iampa';
s.mechanisms(1).equations=CLOCKampa;
data=SimulateModel(s,'time_limits',[0 200],'solver','rk1','dt',.01,'verbose_flag',1);
PlotData(data); set(gcf,'position',[0.3203 0.6241 0.2641 0.1843]);
% Ne=Ni=1: .052 sec

% HH with event-based synapses
s=[];
s.pops(1).name='E';
s.pops(1).size=Ne;
s.pops(1).equations=cat(1,HH,'monitor V.spikes(0)');
s.pops(2).name='I';
s.pops(2).size=Ni;
s.pops(2).equations=cat(1,HH,'monitor V.spikes(0)');
s.pops(2).parameters={'I',0,'noise',0};
s.cons(1).direction='E->I';
s.cons(1).mechanism_list='iampa';
s.mechanisms(1).name='iampa';
s.mechanisms(1).equations=EVENTampa;
data=SimulateModel(s,'time_limits',[0 200],'solver','rk1','dt',.01,'verbose_flag',1);
PlotData(data); set(gcf,'position',[0.3203 0.6241 0.2641 0.1843]);
% Ne=Ni=1: .066 sec

% LIF with clock-based synapses
s=[];
s.pops(1).name='E';
s.pops(1).size=Ne;
s.pops(1).equations=LIF;
s.pops(2).name='I';
s.pops(2).size=Ni;
s.pops(2).equations=LIF;
s.pops(2).parameters={'I',0,'noise',0};
s.cons(1).direction='E->I';
s.cons(1).mechanism_list='iampa';
s.mechanisms(1).name='iampa';
s.mechanisms(1).equations=CLOCKampa;
data=SimulateModel(s,'time_limits',[0 200],'solver','rk1','dt',.01,'verbose_flag',1);
PlotData(data); set(gcf,'position',[0.3203 0.6241 0.2641 0.1843]);
% Ne=Ni=1: .135 sec

% LIF with event-based synapses
s=[];
s.pops(1).name='E';
s.pops(1).size=Ne;
s.pops(1).equations=LIF;
s.pops(2).name='I';
s.pops(2).size=Ni;
s.pops(2).equations=LIF;
s.pops(2).parameters={'I',0,'noise',0};
s.cons(1).direction='E->I';
s.cons(1).mechanism_list='iampa';
s.mechanisms(1).name='iampa';
s.mechanisms(1).equations=EVENTampa;
data=SimulateModel(s,'time_limits',[0 200],'solver','rk1','dt',.01,'verbose_flag',1);
PlotData(data); set(gcf,'position',[0.3203 0.6241 0.2641 0.1843]);
% Ne=Ni=1: .255 sec


%% LIF network: clock-based vs event-based conditional synapses with full connectivity
% i.e., do matrix multiplication (s*netcon) for all elements only when spikes occur
% issue: limiting matrix multiplication to event-related time steps

% iAMPA.mech:
% max_int = 5*tauD; % maximum time over which to look for previous spikes
% if(any(t<tspike_pre+max_int,1))
  % (
  % Isyn(V)=gSYN.*(sum(f(t-tspike_pre-delay))*netcon).*(V-ESYN)
  % )
  % else
  % (
  % Isyn(V)=0
  % )
%   
% may need to initialize Isyn(V)=0 and remove from model.functions

EVENTampa={
  'gSYN=.1; ESYN=0; tauD=2; tauR=0.4; delay=15; max_int=5*tauD'
  'netcon=ones(N_pre,N_post)'
  'f(x) = (exp(-x/tauD)-exp(-x/tauR)).*(x>0)'
  'if(any(t<tspike_pre+max_int,1))(Isyn(X)=gSYN.*(sum(f(t-tspike_pre-delay))*netcon).*(X-ESYN))(Isyn(X)=0)'
  '@isyn += Isyn(X_post)'
};


%% LIF network: clock-based vs event-based conditional synapses with sparse connectivity
% i.e., do matrix multiplication (s*netcon) for relevant elements only when spikes occur
% issue: updating only the relevant matrix elements

netcon=[0 1;1 0;1 0];

% WriteDynaSimSolver.m:
X=netcon;
if numel(find(X==0))/numel(X)>=.5
  [i,j,s]=find(X);
  [m,n]=size(X);
  X=sparse(i,j,s,m,n);
  netcon=X;
end
issparse(netcon)

spre=[.5 .5 .5];
(spre*netcon)

% index=find(t<tspike_pre+max_int)