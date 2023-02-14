function y = dlModelPredictivePFC(Ne, Ni, Nio, NoiseRate)

    fprintf("\n->Initialization of dlPFC Laminar Model: ");
    fprintf("\n-->Excitatory layer neuronal population size = %d (recommended to be multiple of 6), ", Ne);
    fprintf("\n-->Inhibitory layer neuronal population size = %d (recommended to be multiple of 6 and at least 3 times smaller than Excitatory size), ", Ni);
    fprintf("\n-->Input layer (terminal) neuronal population size = %d (arbitrary, better to be close to Inhibitory layer size), ", Nio);
    fprintf("\n-->Overall noise rate = %.4f, ", NoiseRate);

    k1 = 0.07; % Diff. for normal weights (uniform random)
    k2 = 0.04; % Min connectivity weight
    k3 = 0.17; % Diff. for strengthen weights
    k4 = 0.74; % Min. for strengthen weights

    % Connectivity matrices

    % E->I
    Kei = k1*rand(Ne, Ni) + k2; % All-to-all, connectivity from E cells to I cells; mid, sup, deep
    % I->E
    Kie = k1*rand(Ni, Ne) + k2; % All-to-all, connectivity from I cells to E cells; mid, sup, deep

    % E->E
    Kee = k1*rand(Ne, Ne) + k2; % Recurrent E-to-E: mid, sup, deep
    Kii = k1*rand(Ni, Ni) + k2; % Recurrent I-to-I: mid, sup, deep

    kzio = zeros(Nio, Nio); % Null (zero) matrix for disconnections
    KdeepEI = Kie * 4;

    a1 = 1;a2 = ceil(1*Ne/6);
    b1 = ceil(1 + 1*Ne/6);b2 = ceil(2*Ne/6);
    c1 = ceil(1 + 2*Ne/6);c2 = ceil(3*Ne/6);
    cx1_1 = ceil(1 + 3*Ne/6);cx1_2 = ceil(4*Ne/6);
    cx2_1 = ceil(1 + 4*Ne/6);cx2_2 = ceil(5*Ne/6);

    cnx1_1 = 1;cnx1_2 = ceil(1*Ne/6);
    cnx2_1 = ceil(1 + 1*Ne/6);cnx2_2 = ceil(2*Ne/6);
    cnx3_1 = ceil(1 + 2*Ne/6);cnx3_2 = ceil(3*Ne/6);
    cny1_1 = ceil(1 + 3*Ne/6);cny1_2 = ceil(4*Ne/6);
    cny2_1 = ceil(1 + 4*Ne/6);cny2_2 = ceil(5*Ne/6);
    cny3_1 = ceil(1 + 5*Ne/6);cny3_2 = ceil(6*Ne/6);

    KmidEsupE = Kee * 0.3;
    KmidEsupE(a1:a2, [cnx1_1:cnx1_2, cny1_1:cny1_2]) = k3*rand((a2-a1+1), (cnx1_2 - cnx1_1 + cny1_2 - cny1_1 + 2)) + k4; % A -> X1, Y1
    KmidEsupE(b1:b2, [cnx2_1:cnx2_2, cny2_1:cny2_2]) = k3*rand((a2-a1+1), (cnx2_2 - cnx2_1 + cny2_2 - cny2_1 + 2)) + k4; % B -> X2, Y2
    KmidEsupE(c1:c2, [cnx3_1:cnx3_2, cny3_1:cny3_2]) = k3*rand((a2-a1+1), (cnx3_2 - cnx3_1 + cny3_2 - cny3_1 + 2)) + k4; % C -> X3, Y3
    KmidEsupE(cx1_1:cx1_2, cnx1_1:cnx3_2) = k3*rand((a2-a1+1), (cnx3_2 - cnx1_1 + 1)) + k4; % Cx1 -> X1, X2, X3
    KmidEsupE(cx2_1:cx2_2, cny1_1:cny3_2) = k3*rand((a2-a1+1), (cny3_2 - cny1_1 + 1)) + k4; % Cx2 -> Y1, Y2, Y3

    KsupEdeepE = Kee * 0.3;
    KsupEdeepE(cnx1_1:cnx1_2, a1:a2) = k3*rand((cnx1_2 - cnx1_1 + 1), (a2 - a1 + 1)) + k4; % X1 -> O1
    KsupEdeepE(cnx2_1:cnx2_2, b1:b2) = k3*rand((cnx2_2 - cnx2_1 + 1), (b2 - b1 + 1)) + k4; % X2 -> O2
    KsupEdeepE(cnx3_1:cnx3_2, c1:c2) = k3*rand((cnx3_2 - cnx3_1 + 1), (c2 - c1 + 1)) + k4; % X3 -> O3
%     KsupEdeepE(cny1_1:cny1_2, a1:a2) = k3*rand((cny1_2 - cny1_1 + 1), (a2 - a1 + 1)) + k4; % Y1 -> O3
%     KsupEdeepE(cny2_1:cny2_2, b1:b2) = k3*rand((cny2_2 - cny2_1 + 1), (b2 - b1 + 1)) + k4; % Y2 -> O2
%     KsupEdeepE(cny3_1:cny3_2, c1:c2) = k3*rand((cny3_2 - cny3_1 + 1), (c2 - c1 + 1)) + k4; % Y3 -> O1

    KmidEdeepE = Kee * 0.3;
    KmidIdeepE = Kie * 0.3;

    % Time constants
    tauGABA_gamma = 4.07; % ms, decay time constant of inhibition for gamma (50Hz)
    tauGABA_beta = 37.07; % ms, decay time constant of inhibition for beta (25Hz)
    tauAMPA = 4.47; % ms, decay time constant of fast excitation (AMPA)
    %     tauAMPA_beta = 38.4;

    % Maximal synaptic strengths
    gAMPA_ei = .2*(21/Ne); % E->I within layer
    gAMPA_ffee = .2*(21/Ne); % feedforward E->E, mid->sup, sup->deep
    gGABAa_ffie = .2*(21/Ne); % feedforward I->E, mid->deep
    gAMPA_in = .2*(21/Ne);

    gAMPA_ee = 0.11*(21/Ne); % E->E within layer
    gGABAa_ie = 3*(21/Ne); % I->E within layer
    gGABAa_ii = 0.11*(21/Ne); % I->I within layer

    % neuronal dynamics
    eqns = 'dV/dt = (Iapp + @current + noise*randn(1, Npop))/C; Iapp=0; noise=0; C=1; V(0) = -rand(1, Npop)*74;';
    eqns2 = 'dV/dt = (rand(1) + 4.5)*(20*(exp(- (t - t1).^2) - exp(- (t - t2).^2)) + noise*randn(1, Npop))/C; f1=4; t1=10; t2=100; noise=0; C=1; V(0) = -60 - rand(1, Npop)*17;';

    g_poisson = 6.7e-5;

    % cell type
    %     spn_cells = {'spn_iNa','spn_iK','spn_iLeak','spn_iM','spn_iCa','spn_CaBuffer','spn_iKca', 'ctx_iPoisson'};
    ctx_cells = {'iNa','iK', 'ctx_iPoisson'};

    cell_type = ctx_cells; % choose spn_cells and ctx_cells

    % Structures: PING template
    ping=[];

    % E-cells
    ping.populations(1).name = 'E';
    ping.populations(1).size = Ne;
    ping.populations(1).equations = eqns;
    ping.populations(1).mechanism_list = cell_type;
    ping.populations(1).parameters = {'Iapp', 4,'noise', NoiseRate*2, 'g_poisson',g_poisson,'onset_poisson',0,'offset_poisson',0};

    % I-cells
    ping.populations(2).name = 'I';
    ping.populations(2).size = Ni;
    ping.populations(2).equations = eqns;
    ping.populations(2).mechanism_list = cell_type;
    ping.populations(2).parameters = {'Iapp',0,'noise', NoiseRate, 'g_poisson',g_poisson,'onset_poisson',0,'offset_poisson',0};

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
    IOping.populations(1).equations = eqns2;
    IOping.populations(1).mechanism_list = cell_type;
    IOping.populations(1).parameters = {'f1', 1,'noise', 4, 'g_poisson',g_poisson, 't1', 200, 't2', 200};

    % I-cells
    IOping.populations(2).name = 'I';
    IOping.populations(2).size = Nio;
    IOping.populations(2).equations = eqns2;
    IOping.populations(2).mechanism_list = cell_type;
    IOping.populations(2).parameters = {'f1', 1,'noise', 4, 'g_poisson',g_poisson, 't1', 200, 't2', 200};

    % E/I connectivity
    IOping.connections(1).direction = 'E->I';
    IOping.connections(1).mechanism_list = {'iPoisson'};
    IOping.connections(1).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',kzio};

    IOping.connections(2).direction = 'E->E';
    IOping.connections(2).mechanism_list = {'iPoisson'};
    IOping.connections(2).parameters = {'gAMPA',gAMPA_ee,'tauAMPA',tauAMPA,'netcon',kzio};

    IOping.connections(3).direction = 'I->E';
    IOping.connections(3).mechanism_list = {'iPoisson'};
    IOping.connections(3).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',kzio};

    IOping.connections(4).direction = 'I->I';
    IOping.connections(4).mechanism_list = {'iPoisson'};
    IOping.connections(4).parameters = {'gGABAa',gGABAa_ii,'tauGABA',tauGABA_gamma,'netcon',kzio};

    % create independent layers
    sup = dsApplyModifications(ping,{'E','name','supE'; 'I','name','supI'}); % superficial layer (~gamma)
    mid = dsApplyModifications(ping,{'E','name','midE'; 'I','name','midI'}); % middle layer (~gamma)
    deep = dsApplyModifications(ping,{'E','name','deepE'; 'I','name','deepI'}); % deep layer (~beta)
    stimuli1 = dsApplyModifications(IOping,{'E','name','IO_SA1'; 'I','name','IO_SB1'}); % I/O layer (stimuli)
    stimuli2 = dsApplyModifications(IOping,{'E','name','IO_SC1'; 'I','name','IO_SA2'}); % I/O layer (stimuli)
    stimuli3 = dsApplyModifications(IOping,{'E','name','IO_SB2'; 'I','name','IO_SC2'}); % I/O layer (stimuli)
    contex = dsApplyModifications(IOping,{'E','name','IO_Cx1'; 'I','name','IO_Cx2'}); % I/O layer (contex)

    % update deep layer parameters to produce beta rhythm (25Hz)
    deep = dsApplyModifications(deep,{'deepI->deepE','tauGABA',tauGABA_beta});
    deep = dsApplyModifications(deep,{'deepI->deepE','netcon',KdeepEI});

    % create full cortical specification
    s = dsCombineSpecifications(sup, mid, deep, stimuli1, stimuli2, stimuli3, contex);

    % connect the layers and inputs
    fprintf("\n--->Connecting separate layers and inputs:");
    tempconn = zeros(Nio, Ne);

    % Input SA -> midE [1-4]
    Aconn = tempconn;
    Aconn(:, a1:a2) =  0.47;

    c = length(s.connections) + 1;
    s.connections(c).direction = 'IO_SA1->midE';
    s.connections(c).mechanism_list={'iPoisson'};
    s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Aconn};

    c = length(s.connections) + 1;
    s.connections(c).direction = 'IO_SA2->midE';
    s.connections(c).mechanism_list={'iPoisson'};
    s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Aconn};

    % Input SB -> midE [5-8]
    Bconn = tempconn;
    Bconn(:, b1:b2) =  0.47;

    c = length(s.connections)+1;
    s.connections(c).direction = 'IO_SB1->midE';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Bconn};

    c = length(s.connections)+1;
    s.connections(c).direction = 'IO_SB2->midE';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Bconn};

    % Input SC -> midE [9-12]
    Cconn = tempconn;
    Cconn(:, c1:c2) =  0.47;

    c = length(s.connections)+1;
    s.connections(c).direction = 'IO_SC1->midE';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Cconn};

    c = length(s.connections)+1;
    s.connections(c).direction = 'IO_SC2->midE';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Cconn};

    % Contex Cx1 -> midE [13-16]
    Cx1conn = tempconn;
    Cx1conn(:, cx1_1:cx1_2) =  0.47;

    c = length(s.connections)+1;
    s.connections(c).direction = 'IO_Cx1->midE';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Cx1conn};

    % Contex Cx2 -> midE [17-20]
    Cx2conn = tempconn;
    Cx2conn(:, cx2_1:cx2_2) =  0.47;

    c = length(s.connections)+1;
    s.connections(c).direction = 'IO_Cx2->midE';
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

    fprintf("\n->Initialization of dlPFC done. \n");

end