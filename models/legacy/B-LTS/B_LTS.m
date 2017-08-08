% Model: B and LTS from Kramer 2008

NB=1; NLTS=1; Cm=.9; onset=500;
stimB=0;    % -16;
stimLTS=8;  % -40
noiseB=0;   % 30
noiseLTS=0; % 50
gBB=0;      % 15;
gBLTS=0;    % 8
gLTSLTS=0;  % 5
gLTSB=0;

cd /home/jason/models/dnsim/Kramer08/B-LTS;
spec=[];
spec.nodes(1).label = 'B';
spec.nodes(1).multiplicity = NB;
spec.nodes(1).dynamics = {'V''=(current)/Cm'};
spec.nodes(1).mechanisms = {'itonic','randn','iNaF','iKDR','leak'};
spec.nodes(1).parameters = {'Cm',Cm,'V_IC',-65,'IC_noise',0,'g_l',.3,'E_l',-75,...
  'stim',stimB,'onset',onset,'noise',noiseB,'gNaF',200,'gKDR',20};
spec.nodes(2).label = 'LTS';
spec.nodes(2).multiplicity = NLTS;
spec.nodes(2).dynamics = {'V''=(current)/Cm'};
spec.nodes(2).mechanisms = {'itonic','randn','iNaF','iKDR','iAR','leak'};
spec.nodes(2).parameters = {'Cm',Cm,'V_IC',-65,'IC_noise',0,'g_l',1,'E_l',-80,...
  'stim',stimLTS,'onset',onset,'noise',noiseLTS,'gNaF',200,'gKDR',10,'gAR',20};
spec.connections(1,1).label = 'B-B';
spec.connections(1,1).mechanisms = {'iSYN'};
spec.connections(1,1).parameters = {'g_SYN',gBB,'E_SYN',-75,'tauDx',5,'tauRx',.5,'IC_noise',0};
spec.connections(1,2).label = 'B-LTS';
spec.connections(1,2).mechanisms = {'iSYN'};
spec.connections(1,2).parameters = {'g_SYN',gBLTS,'E_SYN',-80,'tauDx',6,'tauRx',.5,'IC_noise',0};
spec.connections(2,2).label = 'LTS-LTS';
spec.connections(2,2).mechanisms = {'iSYN'};
spec.connections(2,2).parameters = {'g_SYN',gLTSLTS,'E_SYN',-80,'tauDx',20,'tauRx',.5,'IC_noise',0};
spec.connections(2,2).label = 'LTS-B';
spec.connections(2,2).mechanisms = {'iSYN'};
spec.connections(2,2).parameters = {'g_SYN',gLTSB,'E_SYN',-80,'tauDx',20,'tauRx',.5,'IC_noise',0};

% process specification and simulate model
data = runsim(spec,'timelimits',[0 1000],'dt',.01,'dsfact',10);
plotv(data,spec,'varlabel','V');
title(sprintf('stim=%g',stimLTS));
% dnsim(spec);

%{
% Parse and compare DNSim models:
[ODEFUN,IC,functions,auxvars,FULLSPEC]=buildmodel2(spec,'verbose',1); % parse DNSim spec structure (returns ODEFUN handle and initial conditions IC for manual simulation and FULLSPEC for automatic simulation)
report=modeldiff(FULLSPEC,FULLSPEC) % compare two DNSim models
% DNSim simulation and plots:
sim_data = biosim(FULLSPEC,'timelimits',[0 100],'dt',.02,'dsfact',10); % wrapper to simulate DNSim models
plotv(sim_data,FULLSPEC,'varlabel','V'); % quickly plot select variables
visualizer(sim_data); % ugly interactive tool hijacked to visualize sim_data
%}
% Manual simulation and plots:
%{
%-----------------------------------------------------------
% Auxiliary variables:
	B_itonic_offset      = inf;
	B_B_iSYN_fanout      = inf;
	B_B_iSYN_UB          = max((10),(10));
	B_B_iSYN_Xpre        = linspace(1,B_B_iSYN_UB,(10))'*ones(1,(10));
	B_B_iSYN_Xpost       = (linspace(1,B_B_iSYN_UB,(10))'*ones(1,(10)))';
	B_B_iSYN_mask        = abs(B_B_iSYN_Xpre-B_B_iSYN_Xpost)<=B_B_iSYN_fanout;
	LTS_itonic_offset    = inf;
	B_LTS_iSYN_fanout    = inf;
	B_LTS_iSYN_UB        = max((10),(10));
	B_LTS_iSYN_Xpre      = linspace(1,B_LTS_iSYN_UB,(10))'*ones(1,(10));
	B_LTS_iSYN_Xpost     = (linspace(1,B_LTS_iSYN_UB,(10))'*ones(1,(10)))';
	B_LTS_iSYN_mask      = abs(B_LTS_iSYN_Xpre-B_LTS_iSYN_Xpost)<=B_LTS_iSYN_fanout;
	LTS_LTS_iSYN_fanout  = inf;
	LTS_LTS_iSYN_UB      = max((10),(10));
	LTS_LTS_iSYN_Xpre    = linspace(1,LTS_LTS_iSYN_UB,(10))'*ones(1,(10));
	LTS_LTS_iSYN_Xpost   = (linspace(1,LTS_LTS_iSYN_UB,(10))'*ones(1,(10)))';
	LTS_LTS_iSYN_mask    = abs(LTS_LTS_iSYN_Xpre-LTS_LTS_iSYN_Xpost)<=LTS_LTS_iSYN_fanout;

% Anonymous functions:
	B_itonic_Itonic      = @(t) (-16)*(t>(0) & t<B_itonic_offset); 
	B_iNaF_hinf          = @(B_V) 1./(1+exp((B_V+(58.3))/(6.7)));  
	B_iNaF_htau          = @(B_V) (0.15) + (1.15)./(1+exp((B_V+(37))/(15)));
	B_iNaF_m0            = @(B_V) 1./(1+exp((-B_V-(38))/10));      
	B_iNaF_aH            = @(B_V) (1./(1+exp((B_V+(58.3))/(6.7)))) ./ ((0.15) + (1.15)./(1+exp((B_V+(37))/(15))));
	B_iNaF_bH            = @(B_V) (1-(1./(1+exp((B_V+(58.3))/(6.7)))))./((0.15) + (1.15)./(1+exp((B_V+(37))/(15))));
	B_iNaF_INaF          = @(B_V,h) (200).*(1./(1+exp((-B_V-(38))/10))).^3.*h.*(B_V-(50));
	B_iKDR_minf          = @(B_V) 1./(1+exp((-B_V-(27))/(11.5)));  
	B_iKDR_mtau          = @(B_V) .25+4.35*exp(-abs(B_V+(10))/(10));
	B_iKDR_aM            = @(B_V) (1./(1+exp((-B_V-(27))/(11.5)))) ./ (.25+4.35*exp(-abs(B_V+(10))/(10)));
	B_iKDR_bM            = @(B_V) (1-(1./(1+exp((-B_V-(27))/(11.5)))))./(.25+4.35*exp(-abs(B_V+(10))/(10)));
	B_iKDR_IKDR          = @(B_V,m) (20).*m.^4.*(B_V-(-100));      
	B_B_iSYN_ISYN        = @(B_V,s) ((15).*(s'*B_B_iSYN_mask)'.*(B_V-(-75)));
	LTS_itonic_Itonic    = @(t) (-40)*(t>(0) & t<LTS_itonic_offset);
	LTS_iNaF_hinf        = @(LTS_V) 1./(1+exp((LTS_V+(58.3))/(6.7)));
	LTS_iNaF_htau        = @(LTS_V) (0.15) + (1.15)./(1+exp((LTS_V+(37))/(15)));
	LTS_iNaF_m0          = @(LTS_V) 1./(1+exp((-LTS_V-(38))/10));  
	LTS_iNaF_aH          = @(LTS_V) (1./(1+exp((LTS_V+(58.3))/(6.7)))) ./ ((0.15) + (1.15)./(1+exp((LTS_V+(37))/(15))));
	LTS_iNaF_bH          = @(LTS_V) (1-(1./(1+exp((LTS_V+(58.3))/(6.7)))))./((0.15) + (1.15)./(1+exp((LTS_V+(37))/(15))));
	LTS_iNaF_INaF        = @(LTS_V,h) (200).*(1./(1+exp((-LTS_V-(38))/10))).^3.*h.*(LTS_V-(50));
	LTS_iKDR_minf        = @(LTS_V) 1./(1+exp((-LTS_V-(27))/(11.5)));
	LTS_iKDR_mtau        = @(LTS_V) .25+4.35*exp(-abs(LTS_V+(10))/(10));
	LTS_iKDR_aM          = @(LTS_V) (1./(1+exp((-LTS_V-(27))/(11.5)))) ./ (.25+4.35*exp(-abs(LTS_V+(10))/(10)));
	LTS_iKDR_bM          = @(LTS_V) (1-(1./(1+exp((-LTS_V-(27))/(11.5)))))./(.25+4.35*exp(-abs(LTS_V+(10))/(10)));
	LTS_iKDR_IKDR        = @(LTS_V,m) (10).*m.^4.*(LTS_V-(-100));  
	LTS_iAR_minf         = @(LTS_V) 1 ./ (1+exp(((-87.5)-LTS_V)/(-5.5)));
	LTS_iAR_mtau         = @(LTS_V) 1./((1).*exp(-14.6-.086*LTS_V)+(1).*exp(-1.87+.07*LTS_V));
	LTS_iAR_aM           = @(LTS_V) (1).*((1 ./ (1+exp(((-87.5)-LTS_V)/(-5.5)))) ./ (1./((1).*exp(-14.6-.086*LTS_V)+(1).*exp(-1.87+.07*LTS_V))));
	LTS_iAR_bM           = @(LTS_V) (1).*((1-(1 ./ (1+exp(((-87.5)-LTS_V)/(-5.5)))))./(1./((1).*exp(-14.6-.086*LTS_V)+(1).*exp(-1.87+.07*LTS_V))));
	LTS_iAR_IAR          = @(LTS_V,m) (50).*m.*(LTS_V-(-35));      
	B_LTS_iSYN_ISYN      = @(LTS_V,s) ((8).*(s'*B_LTS_iSYN_mask)'.*(LTS_V-(-80)));
	LTS_LTS_iSYN_ISYN    = @(LTS_V,s) ((5).*(s'*LTS_LTS_iSYN_mask)'.*(LTS_V-(-80)));

% ODE Handle, ICs, integration, and plotting:
ODEFUN = @(t,X) [(((((-16)*(t>(0) & t<B_itonic_offset)))+(((3).*randn((10),1))+((-((200).*(1./(1+exp((-X(1:10)-(38))/10))).^3.*X(11:20).*(X(1:10)-(50))))+((-((20).*X(21:30).^4.*(X(1:10)-(-100))))+((-(((15).*(X(31:40)'*B_B_iSYN_mask)'.*(X(1:10)-(-75)))))+0))))))/(0.9);((1./(1+exp((X(1:10)+(58.3))/(6.7)))) ./ ((0.15) + (1.15)./(1+exp((X(1:10)+(37))/(15))))).*(1-X(11:20))-((1-(1./(1+exp((X(1:10)+(58.3))/(6.7)))))./((0.15) + (1.15)./(1+exp((X(1:10)+(37))/(15))))).*X(11:20);((1./(1+exp((-X(1:10)-(27))/(11.5)))) ./ (.25+4.35*exp(-abs(X(1:10)+(10))/(10)))).*(1-X(21:30))-((1-(1./(1+exp((-X(1:10)-(27))/(11.5)))))./(.25+4.35*exp(-abs(X(1:10)+(10))/(10)))).*X(21:30);-X(31:40)./(5) + ((1-X(31:40))/(0.5)).*(1+tanh(X(1:10)/10));(((((-40)*(t>(0) & t<LTS_itonic_offset)))+(((5).*randn((10),1))+((-((200).*(1./(1+exp((-X(41:50)-(38))/10))).^3.*X(51:60).*(X(41:50)-(50))))+((-((10).*X(61:70).^4.*(X(41:50)-(-100))))+((-((50).*X(71:80).*(X(41:50)-(-35))))+((-(((8).*(X(81:90)'*B_LTS_iSYN_mask)'.*(X(41:50)-(-80)))))+((-(((5).*(X(91:100)'*LTS_LTS_iSYN_mask)'.*(X(41:50)-(-80)))))+0))))))))/(0.9);((1./(1+exp((X(41:50)+(58.3))/(6.7)))) ./ ((0.15) + (1.15)./(1+exp((X(41:50)+(37))/(15))))).*(1-X(51:60))-((1-(1./(1+exp((X(41:50)+(58.3))/(6.7)))))./((0.15) + (1.15)./(1+exp((X(41:50)+(37))/(15))))).*X(51:60);((1./(1+exp((-X(41:50)-(27))/(11.5)))) ./ (.25+4.35*exp(-abs(X(41:50)+(10))/(10)))).*(1-X(61:70))-((1-(1./(1+exp((-X(41:50)-(27))/(11.5)))))./(.25+4.35*exp(-abs(X(41:50)+(10))/(10)))).*X(61:70);((1).*((1 ./ (1+exp(((-87.5)-X(41:50))/(-5.5)))) ./ (1./((1).*exp(-14.6-.086*X(41:50))+(1).*exp(-1.87+.07*X(41:50)))))).*(1-X(71:80))-((1).*((1-(1 ./ (1+exp(((-87.5)-X(41:50))/(-5.5)))))./(1./((1).*exp(-14.6-.086*X(41:50))+(1).*exp(-1.87+.07*X(41:50)))))).*X(71:80);-X(81:90)./(6) + ((1-X(81:90))/(0.5)).*(1+tanh(X(1:10)/10));-X(91:100)./(20) + ((1-X(91:100))/(0.5)).*(1+tanh(X(41:50)/10));];
IC = [-65          -65          -65          -65          -65          -65          -65          -65          -65          -65          0.5          0.5          0.5          0.5          0.5          0.5          0.5          0.5          0.5          0.5         0.34         0.34         0.34         0.34         0.34         0.34         0.34         0.34         0.34         0.34          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          0.1          -65          -65          -65          -65          -65          -65          -65          -65          -65          -65          0.5          0.5          0.5          0.5          0.5          0.5          0.5          0.5          0.5          0.5         0.34         0.34         0.34         0.34         0.34         0.34         0.34         0.34         0.34         0.34          0.3          0.3          0.3          0.3          0.3          0.3          0.3          0.3          0.3          0.3     0.106847     0.101502     0.102233     0.109718     0.105652      0.10276     0.103969     0.103579     0.100309     0.106956     0.107368     0.109824     0.100833     0.109926     0.101627     0.104962     0.108578     0.106621     0.108355     0.102604];

[t,y]=ode23(ODEFUN,[0 100],IC);   % numerical integration
figure; plot(t,y);           % plot all variables/functions
try legend('B\_V','B\_iNaF\_hNaF','B\_iKDR\_mKDR','B\_B\_iSYN\_sSYNpre','LTS\_V','LTS\_iNaF\_hNaF','LTS\_iKDR\_mKDR','LTS\_iAR\_mAR','B\_LTS\_iSYN\_sSYNpre','LTS\_LTS\_iSYN\_sSYNpre'); end
%-
%}
