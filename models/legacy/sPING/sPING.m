% Model: Sparse Pyramidal-Interneuron-Network-Gamma (PING)
cd(fileparts(mfilename('fullpath')));
clear all
% ---------------------------------------------------------
% Parameters
nE=80; nI=20;
% nE=40; nI=10;
tauI=10; gI=.1; gE=.1; stim=7.5; noise=4;
% ---------------------------------------------------------
spec=[];
spec.nodes(1).label = 'E';
spec.nodes(1).multiplicity = nE;
spec.nodes(1).dynamics = {'V''=(current)./Cm'};
spec.nodes(1).mechanisms = {'K','Na','leak','input','randn'};
spec.nodes(1).parameters = {'Cm',1,'V_IC',-70,'stim',stim,'noise',noise};
spec.nodes(2).label = 'I';
spec.nodes(2).multiplicity = nI;
spec.nodes(2).dynamics = {'V''=(current)./Cm'};
spec.nodes(2).mechanisms = {'K','Na','leak','input','randn'};
spec.nodes(2).parameters = {'Cm',1,'V_IC',-70,'stim',0,'noise',noise};
spec.connections(1,2).label = 'E-I';
spec.connections(1,2).mechanisms = {'AMPA'};
spec.connections(1,2).parameters = {'tauDx',2,'g_SYN',gE*(80/nE)};
spec.connections(2,1).label = 'I-E';
spec.connections(2,1).mechanisms = {'GABAa'};
spec.connections(2,1).parameters = {'tauDx',tauI,'g_SYN',gI*(20/nI)};

% process specification and simulate model
tspan=[0 800]; dt=.05; solver='euler'; % solver='rk2'
data = runsim(spec,'timelimits',tspan,'dt',dt,'dsfact',1,'SOLVER',solver,'coder',1);
plotv(data,spec,'varlabel','V');
plotpow(data,spec,'NFFT',8192); % 16384, 10000
% dnsim(spec);

