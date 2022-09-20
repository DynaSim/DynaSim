function y = dlLaminarCortexNetLWK(ModelParameters, populationName)

    NeSuperficial = ModelParameters.NeSuperficial;
    NSomSuperficial = ModelParameters.NSomSuperficial;
    NPvSuperficial = ModelParameters.NPvSuperficial;

    NeMid = ModelParameters.NeMid;
    NSomMid = 0*ModelParameters.NSomMid; % Kopell model; no SOM in mid
    NPvMid = ModelParameters.NPvMid;

    NeDeep = ModelParameters.NeDeep;
    NSomDeep = ModelParameters.NSomDeep;
    NPvDeep = ModelParameters.NPvDeep;

    Nin = ModelParameters.Nin;
    Nout = ModelParameters.Nout;

    Nstim = ModelParameters.Nstim;
    NoiseRate = ModelParameters.NoiseRate;

    fprintf("\n>Initialization of dlLaminarCortex Model: ");
    fprintf("\n-->As this is a Kopell model, We change/force Number of SOM cells in mid layer to be 0.");

    fprintf("\n->Based on Lee&Whittington&Kopell2013");
    fprintf("\n-->Superficial (L1-3) excitatory neurons count = %d , SOM inhibitory = %d , PV inhibitory = %d ", NeSuperficial, NSomSuperficial, NPvSuperficial); % S
    fprintf("\n-->Middle (L4) excitatory neurons count = %d , SOM inhibitory = %d , PV inhibitory = %d ", NeMid, NSomMid, NPvMid); % M
    fprintf("\n-->Deep (L5-6) excitatory neurons count = %d , SOM inhibitory = %d , PV inhibitory = %d ", NeDeep, NSomDeep, NPvDeep); % D 

    fprintf("\n-->Input connections count (terminal) size = %d ", Nin); % Inputs / Stimuli
    fprintf("\n-->Output connections count (terminal) size = %d ", Nout); % Outputs / Probes
    fprintf("\n-->Overall noise rate = %.4f", NoiseRate); % Randomness / Stochasticity
    fprintf("\n--->Population name is %s", populationName); % Name tag or suffix for all names of this dsModel

    k1 = 0.07; % Diff. for normal weights (uniform random)
    k2 = 0.04; % Min connectivity weight

    NeAvg = (NeSuperficial + NeMid + NeDeep) / 3;
%     NiAvg = (NiSuperficial + NiMid + NiDeep) / 3;

    populationName = ['x', populationName];
    % Connectivity matrices

    % sE->sIsom
    KsupEsupSom = k1*rand(NeSuperficial, NSomSuperficial) + k2;
    % sE->sIpv
    KsupEsupPv = k1*rand(NeSuperficial, NPvSuperficial) + k2;
    % sE->dE
    KsupEdeepE = k1*rand(NeSuperficial, NeDeep) + k2;
    % sIsom->sE
    KsupSomsupE = k1*rand(NSomSuperficial, NeSuperficial) + k2;
    % sIsom->sIpv
    KsupSomsupPv = k1*rand(NSomSuperficial, NPvSuperficial) + k2;

    % sE->sE / sI->sI
    KeeSS = k1*rand(NeSuperficial, NeSuperficial) + k2; % Recurrent E-to-E: mid, sup, deep separately
    KiiSS = k1*rand(NPvSuperficial, NPvSuperficial) + k2; % Recurrent I-to-I: mid, sup, deep separately∂

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
    
    kzio = zeros(Nin, Nin); % Null (zero) matrix for disconnections of Input layer.

    % Sub indices for layer decompostion based in I/O (task specific)
    subLayerIndicesInS = zeros(2, Nin);
    subLayerIndicesInM = zeros(2, Nin);
    subLayerIndicesInD = zeros(2, Nin);

    subLayerIndicesOutS = zeros(2, Nout);
    subLayerIndicesOutM = zeros(2, Nout);
    subLayerIndicesOutD = zeros(2, Nout);

    for i = 1:Nin

        subLayerIndicesInS(1, i) = floor((i-1)*NeSuperficial/Nin) + 1;
        subLayerIndicesInS(2, i) = floor((i)*NeSuperficial/Nin);

        subLayerIndicesInM(1, i) = floor((i-1)*NeMid/Nin) + 1;
        subLayerIndicesInM(2, i) = floor((i)*NeMid/Nin);

        subLayerIndicesInD(1, i) = floor((i-1)*NeDeep/Nin) + 1;
        subLayerIndicesInD(2, i) = floor((i)*NeDeep/Nin);

    end

    for i = 1:Nout

        subLayerIndicesOutS(1, i) = floor((i-1)*NeSuperficial/Nout) + 1;
        subLayerIndicesOutS(2, i) = floor((i)*NeSuperficial/Nout);

        subLayerIndicesOutM(1, i) = floor((i-1)*NeMid/Nout) + 1;
        subLayerIndicesOutM(2, i) = floor((i)*NeMid/Nout);

        subLayerIndicesOutD(1, i) = floor((i-1)*NeDeep/Nout) + 1;
        subLayerIndicesOutD(2, i) = floor((i)*NeDeep/Nout);

    end

    KmidEsupE = 0.3 * (k1*rand(NeMid, NeSuperficial) + k2);    
    KmidEdeepE = 0.3 * (k1*rand(NeMid, NeDeep) + k2);
    KmidIdeepE = 0.3 * (k1*rand(NiMid, NeDeep) + k2);
    KsupEdeepE = 0.3 * (k1*rand(NeSuperficial, NeDeep) + k2);

%     KmidEsupE(subLayerIndicesInM(1, 1):subLayerIndicesInM(2, 1), [subLayerIndicesInS(1, 1):subLayerIndicesInS(2, 1), subLayerIndicesInS(1, 4):subLayerIndicesInS(2, 4)]) = k3*rand((subLayerIndicesInM(2, 1) - subLayerIndicesInM(1, 1) + 1), (subLayerIndicesInS(2, 1) - subLayerIndicesInS(1, 1) + subLayerIndicesInS(2, 4) - subLayerIndicesInS(1, 4) + 2)) + k4; % A -> X1, Y1
%     KmidEsupE(subLayerIndicesInM(1, 2):subLayerIndicesInM(2, 2), [subLayerIndicesInS(1, 2):subLayerIndicesInS(2, 2), subLayerIndicesInS(1, 5):subLayerIndicesInS(2, 5)]) = k3*rand((subLayerIndicesInM(2, 2) - subLayerIndicesInM(1, 2) + 1), (subLayerIndicesInS(2, 2) - subLayerIndicesInS(1, 2) + subLayerIndicesInS(2, 5) - subLayerIndicesInS(1, 5) + 2)) + k4; % B -> X2, Y2
%     KmidEsupE(subLayerIndicesInM(1, 3):subLayerIndicesInM(2, 3), [subLayerIndicesInS(1, 3):subLayerIndicesInS(2, 3), subLayerIndicesInS(1, 6):subLayerIndicesInS(2, 6)]) = k3*rand((subLayerIndicesInM(2, 3) - subLayerIndicesInM(1, 3) + 1), (subLayerIndicesInS(2, 3) - subLayerIndicesInS(1, 3) + subLayerIndicesInS(2, 6) - subLayerIndicesInS(1, 6) + 2)) + k4; % C -> X3, Y3
%     KmidEsupE(subLayerIndicesInM(1, 4):subLayerIndicesInM(2, 4), subLayerIndicesInS(1, 1):subLayerIndicesInS(2, 3)) = k3*rand((subLayerIndicesInM(2, 4) - subLayerIndicesInM(1, 4) + 1), (subLayerIndicesInS(2, 3) - subLayerIndicesInS(1, 1) + 1)) + k4; % Cx1 -> X1, X2, X3
%     KmidEsupE(subLayerIndicesInM(1, 5):subLayerIndicesInM(2, 5), subLayerIndicesInS(1, 4):subLayerIndicesInS(2, 6)) = k3*rand((subLayerIndicesInM(2, 5) - subLayerIndicesInM(1, 5) + 1), (subLayerIndicesInS(2, 6) - subLayerIndicesInS(1, 4) + 1)) + k4; % Cx2 -> Y1, Y2, Y3
% 
%     KsupEdeepE = 0.3 * (k1*rand(NeSuperficial, NeDeep) + k2);
%     KsupEdeepE(subLayerIndicesInM(1, 1):subLayerIndicesInM(2, 1), subLayerIndicesInM(1, 1):subLayerIndicesInM(2, 1)) = k3*rand((subLayerIndicesInM(2, 1) - subLayerIndicesInM(1, 1) + 1), (subLayerIndicesInM(2, 1) - subLayerIndicesInM(1, 1) + 1)) + k4; % X1 -> O1
%     KsupEdeepE(subLayerIndicesInM(1, 2):subLayerIndicesInM(2, 2), subLayerIndicesInM(1, 2):subLayerIndicesInM(2, 2)) = k3*rand((subLayerIndicesInM(2, 2) - subLayerIndicesInM(1, 2) + 1), (subLayerIndicesInM(2, 2) - subLayerIndicesInM(1, 2) + 1)) + k4; % X2 -> O2
%     KsupEdeepE(subLayerIndicesInM(1, 3):subLayerIndicesInM(2, 3), subLayerIndicesInM(1, 3):subLayerIndicesInM(2, 3)) = k3*rand((subLayerIndicesInM(2, 3) - subLayerIndicesInM(1, 3) + 1), (subLayerIndicesInM(2, 3) - subLayerIndicesInM(1, 3) + 1)) + k4; % X3 -> O3
%     KsupEdeepE(cny1_1:cny1_2, a1:a2) = k3*rand((cny1_2 - cny1_1 + 1), (a2 - a1 + 1)) + k4; % Y1 -> O3
%     KsupEdeepE(cny2_1:cny2_2, b1:b2) = k3*rand((cny2_2 - cny2_1 + 1), (b2 - b1 + 1)) + k4; % Y2 -> O2
%     KsupEdeepE(cny3_1:cny3_2, c1:c2) = k3*rand((cny3_2 - cny3_1 + 1), (c2 - c1 + 1)) + k4; % Y3 -> O1

    % Time constants
    tauGABA_gamma = 17.14; % ms, decay time constant of inhibition for gamma (around 50Hz)
    tauGABA_beta = 47.74; % ms, decay time constant of inhibition for beta (around 25Hz)
    tauAMPA = 24.96; % ms, decay time constant of fast excitation (AMPA)
    %     tauAMPA_beta = 38.4;

    % Maximal synaptic strengths
    gAMPA_ei = .2*(21/NeAvg); % E->I within layer
    gAMPA_ffee = .2*(21/NeAvg); % feedforward E->E, mid->sup, sup->deep
    gGABAa_ffie = .2*(21/NeAvg); % feedforward I->E, mid->deep
    gAMPA_in = .2*(21/NeAvg);

    gAMPA_ee = 0.11*(21/NeAvg); % E->E within layer
    gGABAa_ie = 3*(21/NeAvg); % I->E within layer
    gGABAa_ii = 0.11*(21/NeAvg); % I->I within layer

    % neuronal dynamics
    eqns = 'dV/dt = (Iapp + @current + noise*randn(1, Npop))/C; Iapp=0; noise=0; C=1; V(0) = -rand(1, Npop)*74;';
    eqns2 = 'dV/dt = (rand(1) + 4.5)*(20*(exp(- (t - t1).^2) - exp(- (t - t2).^2)) + noise*randn(1, Npop))/C; f1=4; t1=10; t2=100; noise=0; C=1; V(0) = -60 - rand(1, Npop)*17;';

    % cell type
    %     spn_cells = {'spn_iNa','spn_iK','spn_iLeak','spn_iM','spn_iCa','spn_CaBuffer','spn_iKca'};
    ctx_cells = {'iNa','iK'};

    cell_type = ctx_cells; % choose spn_cells and ctx_cells

    % BUILDING
    % Structures: SUPERFICIAL LAYER (1-2)
    pingS=[];

    % E-cells
    pingS.populations(1).name = ['supE', populationName];
    pingS.populations(1).size = NeSuperficial;
    pingS.populations(1).equations = eqns;
    pingS.populations(1).mechanism_list = cell_type;
    pingS.populations(1).parameters = {'Iapp', 4,'noise', NoiseRate*2};

    % I-cells
    pingS.populations(2).name = ['supI', populationName];
    pingS.populations(2).size = NiSuperficial;
    pingS.populations(2).equations = eqns;
    pingS.populations(2).mechanism_list = cell_type;
    pingS.populations(2).parameters = {'Iapp',0,'noise', NoiseRate};

    % E/I connectivity
    pingS.connections(1).direction = [pingS.populations(1).name, '->', pingS.populations(2).name];
    pingS.connections(1).source = pingS.populations(1).name;
    pingS.connections(1).target = pingS.populations(2).name;
    pingS.connections(1).mechanism_list = {'iAMPActx'};
    pingS.connections(1).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',KeiSS};

    pingS.connections(2).direction = [pingS.populations(1).name, '->', pingS.populations(1).name];
    pingS.connections(2).source = pingS.populations(1).name;
    pingS.connections(2).target = pingS.populations(1).name;
    pingS.connections(2).mechanism_list = {'iAMPActx'};
    pingS.connections(2).parameters = {'gAMPA',gAMPA_ee,'tauAMPA',tauAMPA,'netcon',KeeSS};

    pingS.connections(3).direction = [pingS.populations(2).name, '->', pingS.populations(1).name];
    pingS.connections(3).source = pingS.populations(2).name;
    pingS.connections(3).target = pingS.populations(1).name;
    pingS.connections(3).mechanism_list = {'iGABActx'};
    pingS.connections(3).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',KieSS};

    pingS.connections(4).direction = [pingS.populations(2).name, '->', pingS.populations(2).name];
    pingS.connections(4).source = pingS.populations(2).name;
    pingS.connections(4).target = pingS.populations(2).name;
    pingS.connections(4).mechanism_list = {'iGABActx'};
    pingS.connections(4).parameters = {'gGABAa',gGABAa_ii,'tauGABA',tauGABA_gamma,'netcon',KiiSS};

    % Structures: MIDDLE LAYER (3-4)
    pingM=[];

    % E-cells
    pingM.populations(1).name = ['midE', populationName];
    pingM.populations(1).size = NeMid;
    pingM.populations(1).equations = eqns;
    pingM.populations(1).mechanism_list = cell_type;
    pingM.populations(1).parameters = {'Iapp', 4,'noise', NoiseRate*2};

    % I-cells
    pingM.populations(2).name = ['midI', populationName];
    pingM.populations(2).size = NiMid;
    pingM.populations(2).equations = eqns;
    pingM.populations(2).mechanism_list = cell_type;
    pingM.populations(2).parameters = {'Iapp',0,'noise', NoiseRate};

    % E/I connectivity
    pingM.connections(1).direction = [pingM.populations(1).name, '->', pingM.populations(2).name];
    pingM.connections(1).source = pingM.populations(1).name;
    pingM.connections(1).target = pingM.populations(2).name;
    pingM.connections(1).mechanism_list = {'iAMPActx'};
    pingM.connections(1).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',KeiMM};

    pingM.connections(2).direction = [pingM.populations(1).name, '->', pingM.populations(1).name];
    pingM.connections(2).source = pingM.populations(1).name;
    pingM.connections(2).target = pingM.populations(1).name;
    pingM.connections(2).mechanism_list = {'iAMPActx'};
    pingM.connections(2).parameters = {'gAMPA',gAMPA_ee,'tauAMPA',tauAMPA,'netcon',KeeMM};

    pingM.connections(3).direction = [pingM.populations(2).name, '->', pingM.populations(1).name];
    pingM.connections(3).source = pingM.populations(2).name;
    pingM.connections(3).target = pingM.populations(1).name;
    pingM.connections(3).mechanism_list = {'iGABActx'};
    pingM.connections(3).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',KieMM};

    pingM.connections(4).direction = [pingM.populations(2).name, '->', pingM.populations(2).name];
    pingM.connections(4).source = pingM.populations(2).name;
    pingM.connections(4).target = pingM.populations(2).name;
    pingM.connections(4).mechanism_list = {'iGABActx'};
    pingM.connections(4).parameters = {'gGABAa',gGABAa_ii,'tauGABA',tauGABA_gamma,'netcon',KiiMM};
    
    % Structures: DEEP LAYER (5-6)
    pingD=[];

    % E-cells
    pingD.populations(1).name = ['deepE', populationName];
    pingD.populations(1).size = NeDeep;
    pingD.populations(1).equations = eqns;
    pingD.populations(1).mechanism_list = cell_type;
    pingD.populations(1).parameters = {'Iapp', 4,'noise', NoiseRate*2};

    % I-cells
    pingD.populations(2).name = ['deepI', populationName];
    pingD.populations(2).size = NiDeep;
    pingD.populations(2).equations = eqns;
    pingD.populations(2).mechanism_list = cell_type;
    pingD.populations(2).parameters = {'Iapp',0,'noise', NoiseRate};

    % E/I connectivity
    pingD.connections(1).direction = [pingD.populations(1).name, '->', pingD.populations(2).name];
    pingD.connections(1).source = pingD.populations(1).name;
    pingD.connections(1).target = pingD.populations(2).name;
    pingD.connections(1).mechanism_list = {'iAMPActx'};
    pingD.connections(1).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',KeiDD};

    pingD.connections(2).direction = [pingD.populations(1).name, '->', pingD.populations(1).name];
    pingD.connections(2).source = pingD.populations(1).name;
    pingD.connections(2).target = pingD.populations(1).name;
    pingD.connections(2).mechanism_list = {'iAMPActx'};
    pingD.connections(2).parameters = {'gAMPA',gAMPA_ee,'tauAMPA',tauAMPA,'netcon',KeeDD};

    pingD.connections(3).direction = [pingD.populations(2).name, '->', pingD.populations(1).name];
    pingD.connections(3).source = pingD.populations(2).name;
    pingD.connections(3).target = pingD.populations(1).name;
    pingD.connections(3).mechanism_list = {'iGABActx'};
    pingD.connections(3).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_beta,'netcon', KieDD};

    pingD.connections(4).direction = [pingD.populations(2).name, '->', pingD.populations(2).name];
    pingD.connections(4).source = pingD.populations(2).name;
    pingD.connections(4).target = pingD.populations(2).name;
    pingD.connections(4).mechanism_list = {'iGABActx'};
    pingD.connections(4).parameters = {'gGABAa',gGABAa_ii,'tauGABA',tauGABA_beta,'netcon',KiiDD};

    % PING template
    IOping = cell(Nstim, 1);

    for i = 1:Nstim

        % E-cells
        IOping{i}.populations(1).name = ['E', char(64+i), populationName];
        IOping{i}.populations(1).size = Nin;
        IOping{i}.populations(1).equations = eqns2;
        IOping{i}.populations(1).mechanism_list = cell_type;
        IOping{i}.populations(1).parameters = {'f1', 1,'noise', 4};
    
        % I-cells
        IOping{i}.populations(2).name = ['I', char(64+i), populationName];
        IOping{i}.populations(2).size = Nin;
        IOping{i}.populations(2).equations = eqns2;
        IOping{i}.populations(2).mechanism_list = cell_type;
        IOping{i}.populations(2).parameters = {'f1', 1,'noise', 4};
    
        % E/I connectivity
        IOping{i}.connections(1).direction = [IOping{i}.populations(1).name, '->', IOping{i}.populations(2).name];
        IOping{i}.connections(1).source = IOping{i}.populations(1).name;
        IOping{i}.connections(1).target = IOping{i}.populations(2).name;
        IOping{i}.connections(1).mechanism_list = {'iAMPActx'};
        IOping{i}.connections(1).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',kzio};
    
        IOping{i}.connections(2).direction = [IOping{i}.populations(1).name, '->', IOping{i}.populations(1).name];
        IOping{i}.connections(2).source = IOping{i}.populations(1).name;
        IOping{i}.connections(2).target = IOping{i}.populations(1).name;
        IOping{i}.connections(2).mechanism_list = {'iAMPActx'};
        IOping{i}.connections(2).parameters = {'gAMPA',gAMPA_ee,'tauAMPA',tauAMPA,'netcon',kzio};
    
        IOping{i}.connections(3).direction = [IOping{i}.populations(2).name, '->', IOping{i}.populations(1).name];
        IOping{i}.connections(3).source = IOping{i}.populations(2).name;
        IOping{i}.connections(3).target = IOping{i}.populations(1).name;
        IOping{i}.connections(3).mechanism_list = {'iAMPActx'};
        IOping{i}.connections(3).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',kzio};
    
        IOping{i}.connections(4).direction = [IOping{i}.populations(2).name, '->', IOping{i}.populations(2).name];
        IOping{i}.connections(4).source = IOping{i}.populations(2).name;
        IOping{i}.connections(4).target = IOping{i}.populations(2).name;
        IOping{i}.connections(4).mechanism_list = {'iAMPActx'};
        IOping{i}.connections(4).parameters = {'gGABAa',gGABAa_ii,'tauGABA',tauGABA_gamma,'netcon',kzio};
        
    end
    
    % create independent layers
%     sup = dsApplyModifications(pingS,{'E','name','supE'; 'I','name','supI'}); % superficial layer (~gamma)
%     mid = dsApplyModifications(pingM,{'E','name','midE'; 'I','name','midI'}); % middle layer (~gamma)
%     deep = dsApplyModifications(pingD,{'E','name','deepE'; 'I','name','deepI'}); % deep layer (~beta)
%     stimuli1 = dsApplyModifications(IOping,{'E','name','IO_SA1'; 'I','name','IO_SB1'}); % I/O layer (stimuli)
%     stimuli2 = dsApplyModifications(IOping,{'E','name','IO_SC1'; 'I','name','IO_SA2'}); % I/O layer (stimuli)
%     stimuli3 = dsApplyModifications(IOping,{'E','name','IO_SB2'; 'I','name','IO_SC2'}); % I/O layer (stimuli)
%     contex = dsApplyModifications(IOping,{'E','name','IO_Cx1'; 'I','name','IO_Cx2'}); % I/O layer (context)

    % create full cortical specification
    
%     y = stimuli1;
% 
% end

    s = dsCombineSpecifications(pingS, pingM, pingD, IOping{1}, IOping{2}, IOping{3});

    % connect the layers and inputs
    fprintf("\n--->Connecting separate layers and inputs:");
    tempconn = zeros(Nin, NeMid);

    % Inputs/stimuli -> midE

    for k = 1:Nstim

        Aconn = tempconn;
        Aconn(:, subLayerIndicesInM(1, k):subLayerIndicesInM(2, k)) =  0.47;
    
        c = length(s.connections) + 1;
        s.connections(c).direction = [IOping{k}.populations(1).name, '->', pingM.populations(1).name];
        s.connections(c).source = IOping{k}.populations(1).name;
        s.connections(c).target = pingM.populations(1).name;
        s.connections(c).mechanism_list={'iAMPActx'};
        s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Aconn};
    
        c = length(s.connections) + 1;
        s.connections(c).direction = [IOping{k}.populations(2).name, '->', pingM.populations(1).name];
        s.connections(c).source = IOping{k}.populations(2).name;
        s.connections(c).target = pingM.populations(1).name;
        s.connections(c).mechanism_list={'iAMPActx'};
        s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Aconn};

    end

    % midE -> supE
    c = length(s.connections)+1;
    s.connections(c).direction = [pingM.populations(1).name, '->', pingS.populations(1).name];
    s.connections(c).source = pingM.populations(1).name;
    s.connections(c).target = pingS.populations(1).name;
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_ffee,'tauAMPA',tauAMPA,'netcon',KmidEsupE};

    % midE -> deepE
    c = length(s.connections)+1;
    s.connections(c).direction = [pingM.populations(1).name, '->', pingD.populations(1).name];
    s.connections(c).source = pingM.populations(1).name;
    s.connections(c).target = pingD.populations(1).name;
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_ffee,'tauAMPA',tauAMPA,'netcon',KmidEdeepE};

    % midI -> deepE
    c = length(s.connections)+1;
    s.connections(c).direction = [pingM.populations(2).name, '->', pingD.populations(1).name];
    s.connections(c).source = pingM.populations(2).name;
    s.connections(c).target = pingD.populations(1).name;
    s.connections(c).mechanism_list={'iGABActx'};
    s.connections(c).parameters={'gGABAa',gGABAa_ffie,'tauGABA',tauGABA_beta,'netcon',KmidIdeepE};

    % supE -> deepE
    c = length(s.connections)+1;
    s.connections(c).direction = [pingS.populations(1).name, '->', pingD.populations(1).name];
    s.connections(c).source = pingS.populations(1).name;
    s.connections(c).target = pingD.populations(1).name;
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_ffee,'tauAMPA',tauAMPA,'netcon',KsupEdeepE};

    % Outputs: deepE is the recommended output layer.
    y = s;
    fprintf("\n->Initialization of dsModel ""%s"" is done. \n", populationName);

end