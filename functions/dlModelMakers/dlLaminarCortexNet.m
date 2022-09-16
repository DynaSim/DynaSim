function y = dlLaminarCortexNet(NeSuperficial, NiSuperficial, NeMid, NiMid, NeDeep, NiDeep, Nio, NoiseRate)

    fprintf("\n->Initialization of dlLaminarCortex Model: ");
    fprintf("\n-->Superficial (L1-2) excitatory neurons count = %d , inhibitory = %d ", NeSuperficial, NiSuperficial); % S
    fprintf("\n-->Middle (L3-4) excitatory neurons count = %d , inhibitory = %d ", NeMid, NiMid); % M
    fprintf("\n-->Middle (L3-4) excitatory neurons count = %d , inhibitory = %d ", NeDeep, NiDeep); % D 
    fprintf("\n-->Input connections count (terminal) size = %d ", Nio); % D
    fprintf("\n-->Overall noise rate = %.4f ", NoiseRate);

    k1 = 0.07; % Diff. for normal weights (uniform random)
    k2 = 0.04; % Min connectivity weight
    k3 = 0.17; % Diff. for strengthen weights
    k4 = 0.74; % Min. for strengthen weights

    % Connectivity matrices

    % sE->sI
    KeiSS = k1*rand(NeSuperficial, NiSuperficial) + k2; % All-to-all, connectivity from E cells to I cells; mid, sup, deep
    % sI->sE
    KieSS = k1*rand(NiSuperficial, NeSuperficial) + k2; % All-to-all, connectivity from I cells to E cells; mid, sup, deep

    % sE->sE / sI->sI
    KeeSS = k1*rand(NeSuperficial, NeSuperficial) + k2; % Recurrent E-to-E: mid, sup, deep separately
    KiiSS = k1*rand(NiSuperficial, NiSuperficial) + k2; % Recurrent I-to-I: mid, sup, deep separately

    % mE->mI
    KeiMM = k1*rand(NeMid, NiMid) + k2; % All-to-all, connectivity from E cells to I cells; mid, sup, deep
    % mI->mE
    KieMM = k1*rand(NiMid, NeMid) + k2; % All-to-all, connectivity from I cells to E cells; mid, sup, deep

    % mE->mE / mI->mI
    KeeMM = k1*rand(NeMid, NeMid) + k2; % Recurrent E-to-E: mid, sup, deep separately
    KiiMM = k1*rand(NiMid, NiMid) + k2; % Recurrent I-to-I: mid, sup, deep separately

    % dE->dI
    KeiDD = k1*rand(NeDeep, NiDeep) + k2; % All-to-all, connectivity from E cells to I cells; mid, sup, deep
    % dI->dE
    KieDD = 4*k1*rand(NiDeep, NeDeep) + k2; % All-to-all, connectivity from I cells to E cells; mid, sup, deep

    % dE->dE / dI->dI
    KeeDD = k1*rand(NeDeep, NeDeep) + k2; % Recurrent E-to-E: mid, sup, deep separately
    KiiDD = k1*rand(NiDeep, NiDeep) + k2; % Recurrent I-to-I: mid, sup, deep separately
    
    kzio = zeros(Nio, Nio); % Null (zero) matrix for disconnections

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

    % BUILDING
    % Structures: SUPERFICIAL LAYER (1-2)
    pingS=[];

    % E-cells
    pingS.populations(1).name = 'E';
    pingS.populations(1).size = NeSuperficial;
    pingS.populations(1).equations = eqns;
    pingS.populations(1).mechanism_list = cell_type;
    pingS.populations(1).parameters = {'Iapp', 4,'noise', NoiseRate*2, 'g_poisson',g_poisson,'onset_poisson',0,'offset_poisson',0};

    % I-cells
    pingS.populations(2).name = 'I';
    pingS.populations(2).size = NiSuperficial;
    pingS.populations(2).equations = eqns;
    pingS.populations(2).mechanism_list = cell_type;
    pingS.populations(2).parameters = {'Iapp',0,'noise', NoiseRate, 'g_poisson',g_poisson,'onset_poisson',0,'offset_poisson',0};

    % E/I connectivity
    pingS.connections(1).direction = 'E->I';
    pingS.connections(1).mechanism_list = {'iAMPActx'};
    pingS.connections(1).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',KeiSS};

    pingS.connections(2).direction = 'E->E';
    pingS.connections(2).mechanism_list = {'iAMPActx'};
    pingS.connections(2).parameters = {'gAMPA',gAMPA_ee,'tauAMPA',tauAMPA,'netcon',KeeSS};

    pingS.connections(3).direction = 'I->E';
    pingS.connections(3).mechanism_list = {'iGABActx'};
    pingS.connections(3).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',KieSS};

    pingS.connections(4).direction = 'I->I';
    pingS.connections(4).mechanism_list = {'iGABActx'};
    pingS.connections(4).parameters = {'gGABAa',gGABAa_ii,'tauGABA',tauGABA_gamma,'netcon',KiiSS};

    % Structures: MIDDLE LAYER (3-4)
    pingM=[];

    % E-cells
    pingM.populations(1).name = 'E';
    pingM.populations(1).size = NeMid;
    pingM.populations(1).equations = eqns;
    pingM.populations(1).mechanism_list = cell_type;
    pingM.populations(1).parameters = {'Iapp', 4,'noise', NoiseRate*2, 'g_poisson',g_poisson,'onset_poisson',0,'offset_poisson',0};

    % I-cells
    pingM.populations(2).name = 'I';
    pingM.populations(2).size = NiMid;
    pingM.populations(2).equations = eqns;
    pingM.populations(2).mechanism_list = cell_type;
    pingM.populations(2).parameters = {'Iapp',0,'noise', NoiseRate, 'g_poisson',g_poisson,'onset_poisson',0,'offset_poisson',0};

    % E/I connectivity
    pingM.connections(1).direction = 'E->I';
    pingM.connections(1).mechanism_list = {'iAMPActx'};
    pingM.connections(1).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',KeiMM};

    pingM.connections(2).direction = 'E->E';
    pingM.connections(2).mechanism_list = {'iAMPActx'};
    pingM.connections(2).parameters = {'gAMPA',gAMPA_ee,'tauAMPA',tauAMPA,'netcon',KeeMM};

    pingM.connections(3).direction = 'I->E';
    pingM.connections(3).mechanism_list = {'iGABActx'};
    pingM.connections(3).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',KieMM};

    pingM.connections(4).direction = 'I->I';
    pingM.connections(4).mechanism_list = {'iGABActx'};
    pingM.connections(4).parameters = {'gGABAa',gGABAa_ii,'tauGABA',tauGABA_gamma,'netcon',KiiMM};
    
    % Structures: DEEP LAYER (5-6)
    pingD=[];

    % E-cells
    pingD.populations(1).name = 'E';
    pingD.populations(1).size = NeDeep;
    pingD.populations(1).equations = eqns;
    pingD.populations(1).mechanism_list = cell_type;
    pingD.populations(1).parameters = {'Iapp', 4,'noise', NoiseRate*2, 'g_poisson',g_poisson,'onset_poisson',0,'offset_poisson',0};

    % I-cells
    pingD.populations(2).name = 'I';
    pingD.populations(2).size = NiDeep;
    pingD.populations(2).equations = eqns;
    pingD.populations(2).mechanism_list = cell_type;
    pingD.populations(2).parameters = {'Iapp',0,'noise', NoiseRate, 'g_poisson',g_poisson,'onset_poisson',0,'offset_poisson',0};

    % E/I connectivity
    pingD.connections(1).direction = 'E->I';
    pingD.connections(1).mechanism_list = {'iAMPActx'};
    pingD.connections(1).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',KeiDD};

    pingD.connections(2).direction = 'E->E';
    pingD.connections(2).mechanism_list = {'iAMPActx'};
    pingD.connections(2).parameters = {'gAMPA',gAMPA_ee,'tauAMPA',tauAMPA,'netcon',KeeDD};

    pingD.connections(3).direction = 'I->E';
    pingD.connections(3).mechanism_list = {'iGABActx'};
    pingD.connections(3).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon', KieDD};

    pingD.connections(4).direction = 'I->I';
    pingD.connections(4).mechanism_list = {'iGABActx'};
    pingD.connections(4).parameters = {'gGABAa',gGABAa_ii,'tauGABA',tauGABA_gamma,'netcon',KiiDD};

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
    sup = dsApplyModifications(pingS,{'E','name','supE'; 'I','name','supI'}); % superficial layer (~gamma)
    mid = dsApplyModifications(pingM,{'E','name','midE'; 'I','name','midI'}); % middle layer (~gamma)
    deep = dsApplyModifications(pingD,{'E','name','deepE'; 'I','name','deepI'}); % deep layer (~beta)
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