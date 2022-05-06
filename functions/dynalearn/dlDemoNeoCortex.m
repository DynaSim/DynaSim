function y = dlDemoNeoCortex(Ne, Ni, Nio, noise_rate)

    fprintf("Initialization...\n");

    % Population sizes

%     Ne = 20;     % # of E cells per layer
%     Ni = 4;  % # of I cells per layer
%     Nio = 10; % # of Input cells

    k1 = 0.1; % Diff. for normal weights (uniform random)
    k2 = 0.2; % Min connectivity weight
    k3 = 0.1; % Diff. for strengthen weights
    k4 = 0.9; % Min. for strengthen weights

    % Connectivity matrices

    % E->I
    Kei = k1*rand(Ne, Ni) + 0.3; % all-to-all, connectivity from E cells to I cells; mid, sup, deep
    % I->E
    Kie = k1*rand(Ni, Ne) + 0.3; % all-to-all, connectivity from I cells to E cells; mid, sup, deep

    % E->E
    Kee = k1*rand(Ne, Ne) + k2; % recurrent E-to-E: mid, sup, deep
    Kii = k1*rand(Ni, Ni) + k2; % recurrent I-to-I: mid, sup, deep
%     Kffee = k1*rand(Ne, Ne) + k2; % feedforward E-to-E: mid->sup, sup->deep
%     Kffie = k1*rand(Ni, Ne) + k2; % feedforward I-to-E: mid->deep

    kzio = zeros(Nio, Nio);
    KdeepEI = Kie * 1.5;

    a1 = 1;a2 = ceil(1*Ne/4);
    b1 = ceil(1 + 1*Ne/4);b2 = ceil(2*Ne/4);
    c1 = ceil(1 + 2*Ne/4);c2 = ceil(3*Ne/4);
    d1 = ceil(1 + 3*Ne/4);d2 = ceil(4*Ne/4);
    
    % Manual weight adjustment
%     KmidEmidI = Kei * 0.3;
    % KmidEmidI(1:5, [1, 2]) = k3*rand(5, 2) + k4; % !A -> Z1, Z2
    % KmidEmidI(6:10, [3, 4]) = k3*rand(5, 2) + k4; % !B -> Z3, Z4
    % KmidEmidI(11:15, [1, 4]) = k3*rand(5, 2) + k4; % !C1 -> Z1, Z4
    % KmidEmidI(16:20, [2, 3]) = k3*rand(5, 2) + k4; % !C2 -> Z2, Z3
    
    KmidEsupE = Kee * 0.3;
    KmidEsupE(a1:a2, [a1:a2, c1:c2]) = k3*rand((a2-a1+1), (b2-a1+1)) + k4; % A -> X1, Y1
    KmidEsupE(b1:b2, [b1:b2, d1:d2]) = k3*rand((a2-a1+1), (b2-a1+1)) + k4; % B -> X2, Y2
    KmidEsupE(c1:c2, a1:b2) = k3*rand((a2-a1+1), (b2-a1+1)) + k4; % C1 -> X1, X2
    KmidEsupE(d1:d2, c1:d2) = k3*rand((a2-a1+1), (b2-a1+1)) + k4; % C2 -> Y1, Y2

    KmidEdeepE = Kee * 0.3;
    % KmidEdeepE(1:5, 1:10) = k3*rand((a2-a1+1), (b2-a1+1)) + k4; % A -> O1
    % KmidEdeepE(6:10, 11:20) = k3*rand((a2-a1+1), (b2-a1+1)) + k4; % B -> O2

    KsupEdeepE = Kee*0.3;
    KsupEdeepE(a1:a2, a1:b2) = k3*rand((a2-a1+1), (b2-a1+1)) + k4; % X1 -> O1
    KsupEdeepE(b1:b2, c1:d2) = k3*rand((a2-a1+1), (b2-a1+1)) + k4; % X2 -> O2
    KsupEdeepE(c1:c2, c1:d2) = k3*rand((a2-a1+1), (b2-a1+1)) + k4; % Y1 -> O2
    KsupEdeepE(d1:d2, a1:b2) = k3*rand((a2-a1+1), (b2-a1+1)) + k4; % Y2 -> O1

    KmidIdeepE = Kie*0.3;
    KmidIdeepE(1, c1:d2) = 0.1*rand(1, (b2-a1+1)) + 0.6; % !(A & C1) -> O2 
    KmidIdeepE(2, a1:b2) = 0.1*rand(1, (b2-a1+1)) + 0.6; % !(A & C2) -> O1
    KmidIdeepE(3, c1:d2) = 0.1*rand(1, (b2-a1+1)) + 0.6; % !(B & C2) -> O2
    KmidIdeepE(4, a1:b2) = 0.1*rand(1, (b2-a1+1)) + 0.6; % !(B & C1) -> O1

    % Time constants
    tauGABA_gamma = 4.8; % ms, decay time constant of inhibition for gamma (50Hz)
    tauGABA_beta = 38.4; % ms, decay time constant of inhibition for beta (25Hz)
    tauAMPA = 4.8; % ms, decay time constant of fast excitation (AMPA)
%     tauAMPA_beta = 38.4;

    % Maximal synaptic strengths
    gAMPA_ei = .2*(20/Ne); % E->I within layer
    gAMPA_ffee = .2*(20/Ne); % feedforward E->E, mid->sup, sup->deep
    gGABAa_ffie = .2*(20/Ne); % feedforward I->E, mid->deep
    gAMPA_in = .2*(20/Ne);

    gAMPA_ee = 0.11*(20/Ne); % E->E within layer
    gGABAa_ie = 4*(20/Ne); % I->E within layer
    gGABAa_ii = 0.11*(20/Ne); % I->I within layer
%     noise_rate = 12;

    % neuronal dynamics
    eqns = 'dV/dt = (Iapp + @current + noise*randn(1,Npop))/C; Iapp=0; noise=0; C=1';

    % SPN
%     g_l_D1 = 0.096;      % mS/cm^2, Leak conductance for D1 SPNs 
%     g_l_D2 = 0.1;        % mS/cm^2, Leak conductance for D2 SPNs
%     g_cat_D1 = 0.018;    % mS/cm^2, Conductance of the T-type Ca2+ current for D1 SPNs
%     g_cat_D2 = 0.025;    % mS/cm^2, Conductance of the T-type Ca2+ current for D2 SPNs

    g_poisson = 6.4e-4;

    % cell type
%     spn_cells = {'spn_iNa','spn_iK','spn_iLeak','spn_iM','spn_iCa','spn_CaBuffer','spn_iKca', 'ctx_iPoisson'};
    ctx_cells = {'iNa','iK', 'ctx_iPoisson'};

    cell_type = ctx_cells; % choose spn_cells and ctx_cells

    % create DynaSim specification structure

    % PING template
    ping=[];

    % E-cells
    ping.populations(1).name = 'E';
    ping.populations(1).size = Ne;
    ping.populations(1).equations = eqns;
    ping.populations(1).mechanism_list = cell_type;
    ping.populations(1).parameters = {'Iapp', 3,'noise', noise_rate*3, 'g_poisson',g_poisson,'onset_poisson',0,'offset_poisson',0};

    % I-cells
    ping.populations(2).name = 'I';
    ping.populations(2).size = Ni;
    ping.populations(2).equations = eqns;
    ping.populations(2).mechanism_list = cell_type;
    ping.populations(2).parameters = {'Iapp',0,'noise', noise_rate, 'g_poisson',g_poisson,'onset_poisson',0,'offset_poisson',0};

    % E/I connectivity
    ping.connections(1).direction = 'E->I';
    ping.connections(1).mechanism_list = {'iAMPActx'};
    ping.connections(1).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',Kei};

    ping.connections(2).direction = 'E->E';
    ping.connections(2).mechanism_list = {'iAMPActx'};
    ping.connections(2).parameters = {'gAMPA',gAMPA_ee,'tauAMPA',tauAMPA,'netcon',Kee};

    ping.connections(3).direction = 'I->E';
    ping.connections(3).mechanism_list = {'iGABActx'};
    ping.connections(3).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',Kie};

    ping.connections(4).direction = 'I->I';
    ping.connections(4).mechanism_list = {'iGABActx'};
    ping.connections(4).parameters = {'gGABAa',gGABAa_ii,'tauGABA',tauGABA_gamma,'netcon',Kii};

    % PING template
    IOping=[];

    % E-cells
    IOping.populations(1).name = 'E';
    IOping.populations(1).size = Nio;
    IOping.populations(1).equations = eqns;
    IOping.populations(1).mechanism_list = cell_type;
    IOping.populations(1).parameters = {'Iapp',1,'noise', 2, 'g_poisson',g_poisson,'onset_poisson',0,'offset_poisson',0};

    % I-cells
    IOping.populations(2).name = 'I';
    IOping.populations(2).size = Nio;
    IOping.populations(2).equations = eqns;
    IOping.populations(2).mechanism_list = cell_type;
    IOping.populations(2).parameters = {'Iapp',1,'noise', 2, 'g_poisson',g_poisson,'onset_poisson',0,'offset_poisson',0};

    % E/I connectivity
    IOping.connections(1).direction = 'E->I';
    IOping.connections(1).mechanism_list = {'iAMPActx'};
    IOping.connections(1).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',kzio};

    IOping.connections(2).direction = 'E->E';
    IOping.connections(2).mechanism_list = {'iAMPActx'};
    IOping.connections(2).parameters = {'gAMPA',gAMPA_ee,'tauAMPA',tauAMPA,'netcon',kzio};

    IOping.connections(3).direction = 'I->E';
    IOping.connections(3).mechanism_list = {'iGABActx'};
    IOping.connections(3).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',kzio};

    IOping.connections(4).direction = 'I->I';
    IOping.connections(4).mechanism_list = {'iGABActx'};
    IOping.connections(4).parameters = {'gGABAa',gGABAa_ii,'tauGABA',tauGABA_gamma,'netcon',kzio};

    % create independent layers
    sup = dsApplyModifications(ping,{'E','name','supE'; 'I','name','supI'}); % superficial layer (~gamma)
    mid = dsApplyModifications(ping,{'E','name','midE'; 'I','name','midI'}); % middle layer (~gamma)
    deep = dsApplyModifications(ping,{'E','name','deepE'; 'I','name','deepI'}); % deep layer (~beta)
    stimuli = dsApplyModifications(IOping,{'E','name','SA'; 'I','name','SB'}); % I/O layer (stimuli)
    contex = dsApplyModifications(IOping,{'E','name','Cx1'; 'I','name','Cx2'}); % I/O layer (contex)

    % uppdate deep layer parameters to produce beta rhythm (25Hz)
    deep = dsApplyModifications(deep,{'deepI->deepE','tauGABA',tauGABA_beta});
    deep = dsApplyModifications(deep,{'deepI->deepE','netcon',KdeepEI});

    % create full cortical specification
    s = dsCombineSpecifications(sup, mid, deep, stimuli, contex);

    % connect the layers and inputs

    fprintf("Connecting separate layers and inputs...\n");

    % Input SA -> midE [1-3]
    tempconn = zeros(Nio, Ne);
    Aconn = tempconn;
    Aconn(:, a1:a2) =  1;

    c = length(s.connections) + 1;
    s.connections(c).direction = 'SA->midE';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Aconn};

    % Input SB -> midE [4-6]
    Bconn = tempconn;
    Bconn(:, b1:b2) =  1;
    c = length(s.connections)+1;
    s.connections(c).direction = 'SB->midE';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Bconn};

    % Contex Cx1 -> midE [7-9]
    Cx1conn = tempconn;
    Cx1conn(:, c1:c2) =  1;

    c = length(s.connections)+1;
    s.connections(c).direction = 'Cx1->midE';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Cx1conn};

    % Contex Cx2 -> midE [10-12]
    Cx2conn = tempconn;
    Cx2conn(:, d1:d2) =  1;

    c = length(s.connections)+1;
    s.connections(c).direction = 'Cx2->midE';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Cx2conn};

    % midE -> supE
    c = length(s.connections)+1;
    s.connections(c).direction = 'midE->supE';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_ffee,'tauAMPA',tauAMPA,'netcon',KmidEsupE};
    
    % midE -> deepE
    c = length(s.connections)+1;
    s.connections(c).direction = 'midE->deepE';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_ffee,'tauAMPA',tauAMPA,'netcon',KmidEdeepE};

    % midI -> deepE
    c = length(s.connections)+1;
    s.connections(c).direction = 'midI->deepE';
    s.connections(c).mechanism_list={'iGABActx'};
    s.connections(c).parameters={'gGABAa',gGABAa_ffie,'tauGABA',tauGABA_beta,'netcon',KmidIdeepE};

    % supE -> deepE
    c = length(s.connections)+1;
    s.connections(c).direction = 'supE->deepE';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_ffee,'tauAMPA',tauAMPA,'netcon',KsupEdeepE};

    % Outputs: deepE [1-Ne/2] as O1
    % deepE [Ne/2+1-Ne] as O2
    y = s;
    
    fprintf("Initialization done.\n");

end
