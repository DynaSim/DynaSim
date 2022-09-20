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

    % sIpv->sE
    KsupPvsupE = k1*rand(NPvSuperficial, NeSuperficial) + k2;
    % sIpv->sIsom
    KsupPvsupSom = k1*rand(NPvSuperficial, NSomSuperficial) + k2;

    % mE->sIsom
    KmidEsupSom = k1*rand(NeMid, NSomSuperficial) + k2;
    % mE->mIpv
    KmidEmidPv = k1*rand(NeMid, NPvMid) + k2;
    % mE->sE
    KmidEsupE = k1*rand(NeMid, NeSuperficial) + k2;
    % mE->dE
    KmidEdeepE = k1*rand(NeMid, NeDeep) + k2;
    % mE->dIpv
    KmidEdeepPv = k1*rand(NeMid, NPvDeep) + k2;

    % mIpv->sE
    KmidPvsupE = k1*rand(NPvMid, NeSuperficial) + k2;
    % mIpv->mE
    KmidPvmidE = k1*rand(NPvMid, NeMid) + k2;

    % dE->dIsom
    KdeepEdeepSom = k1*rand(NeDeep, NSomDeep) + k2;
    % dE->sIsom
    KdeepEsupSom = k1*rand(NeDeep, NSomSuperficial) + k2;
    % dE->dIpv
    KdeepEdeepPv = k1*rand(NeDeep, NPvDeep) + k2;
    % dE->sIpv
    KdeepEsupPv = k1*rand(NeDeep, NPvSuperficial) + k2;

    % dIpv->dE
    KdeepPvdeepE = k1*rand(NPvDeep, NeDeep) + k2;
    % dIpv->dISom
    KdeepPvdeepSom = k1*rand(NPvDeep, NSomDeep) + k2;

    % dIsom->dE
    KdeepSomdeepE = k1*rand(NSomDeep, NeDeep) + k2;
    % dIsom->dIpv
    KdeepSomsupPv = k1*rand(NSomDeep, NPvDeep) + k2;
    % dIsom->mIpv
    KdeepSommidPv = k1*rand(NSomDeep, NPvMid) + k2;
    
    KnullIO = zeros(Nin, Nin); % Null (zero) matrix for disconnections of Input layer.

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

    % Time constants
    tauGABA_gamma = 17.14; % ms, decay time constant of inhibition for gamma (around 50Hz)
    tauGABA_beta = 47.74; % ms, decay time constant of inhibition for beta (around 25Hz)
    tauAMPA = 24.96; % ms, decay time constant of fast excitation (AMPA)

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
    ctx_cells = {'iNa','iK'};
    cell_type = ctx_cells;

    % BUILDING
    % Structures: SUPERFICIAL LAYER (1-2-3)
    pingS=[];

    % E-cells
    pingS.populations(1).name = ['supE', populationName];
    pingS.populations(1).size = NeSuperficial;
    pingS.populations(1).equations = eqns;
    pingS.populations(1).mechanism_list = cell_type;
    pingS.populations(1).parameters = {'Iapp', 4,'noise', NoiseRate*2};

    % I-cells-SOM
    pingS.populations(2).name = ['supISOM', populationName];
    pingS.populations(2).size = NSomSuperficial;
    pingS.populations(2).equations = eqns;
    pingS.populations(2).mechanism_list = cell_type;
    pingS.populations(2).parameters = {'Iapp',0,'noise', NoiseRate};

    % I-cells-PV
    pingS.populations(3).name = ['supIPV', populationName];
    pingS.populations(3).size = NPvSuperficial;
    pingS.populations(3).equations = eqns;
    pingS.populations(3).mechanism_list = cell_type;
    pingS.populations(3).parameters = {'Iapp',0,'noise', NoiseRate};
    
    % Interlayer connections - Superficial
    pingS.connections(1).direction = [pingS.populations(1).name, '->', pingS.populations(2).name];
    pingS.connections(1).source = pingS.populations(1).name;
    pingS.connections(1).target = pingS.populations(2).name;
    pingS.connections(1).mechanism_list = {'iAMPActx'};
    pingS.connections(1).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',KsupEsupSom};

    pingS.connections(2).direction = [pingS.populations(1).name, '->', pingS.populations(3).name];
    pingS.connections(2).source = pingS.populations(1).name;
    pingS.connections(2).target = pingS.populations(3).name;
    pingS.connections(2).mechanism_list = {'iAMPActx'};
    pingS.connections(2).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',KsupEsupPv};

    pingS.connections(3).direction = [pingS.populations(1).name, '->', pingS.populations(2).name];
    pingS.connections(3).source = pingS.populations(1).name;
    pingS.connections(3).target = pingS.populations(2).name;
    pingS.connections(3).mechanism_list = {'iNMDActx'};
    pingS.connections(3).parameters = {'netcon',KsupEsupSom};

    pingS.connections(4).direction = [pingS.populations(1).name, '->', pingS.populations(3).name];
    pingS.connections(4).source = pingS.populations(1).name;
    pingS.connections(4).target = pingS.populations(3).name;
    pingS.connections(4).mechanism_list = {'iNMDActx'};
    pingS.connections(4).parameters = {'netcon',KsupEsupPv};

    pingS.connections(5).direction = [pingS.populations(2).name, '->', pingS.populations(1).name];
    pingS.connections(5).source = pingS.populations(2).name;
    pingS.connections(5).target = pingS.populations(1).name;
    pingS.connections(5).mechanism_list = {'iGABActx'};
    pingS.connections(5).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',KsupSomsupE};

    pingS.connections(6).direction = [pingS.populations(2).name, '->', pingS.populations(3).name];
    pingS.connections(6).source = pingS.populations(2).name;
    pingS.connections(6).target = pingS.populations(3).name;
    pingS.connections(6).mechanism_list = {'iGABActx'};
    pingS.connections(6).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',KsupSomsupPv};

    pingS.connections(7).direction = [pingS.populations(3).name, '->', pingS.populations(1).name];
    pingS.connections(7).source = pingS.populations(3).name;
    pingS.connections(7).target = pingS.populations(1).name;
    pingS.connections(7).mechanism_list = {'iGABActx'};
    pingS.connections(7).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',KsupPvsupE};

    pingS.connections(8).direction = [pingS.populations(3).name, '->', pingS.populations(2).name];
    pingS.connections(8).source = pingS.populations(3).name;
    pingS.connections(8).target = pingS.populations(2).name;
    pingS.connections(8).mechanism_list = {'iGABActx'};
    pingS.connections(8).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',KsupPvsupSom};

    % Structures: MIDDLE LAYER (4)
    pingM=[];

    % E-cells
    pingM.populations(1).name = ['midE', populationName];
    pingM.populations(1).size = NeMid;
    pingM.populations(1).equations = eqns;
    pingM.populations(1).mechanism_list = cell_type;
    pingM.populations(1).parameters = {'Iapp', 4,'noise', NoiseRate*2};

    % I-cells-PV
    pingM.populations(2).name = ['midIPV', populationName];
    pingM.populations(2).size = NPvMid;
    pingM.populations(2).equations = eqns;
    pingM.populations(2).mechanism_list = cell_type;
    pingM.populations(2).parameters = {'Iapp',0,'noise', NoiseRate};

    % Interlayer connections - Middle
    pingM.connections(1).direction = [pingM.populations(1).name, '->', pingM.populations(2).name];
    pingM.connections(1).source = pingM.populations(1).name;
    pingM.connections(1).target = pingM.populations(2).name;
    pingM.connections(1).mechanism_list = {'iAMPActx'};
    pingM.connections(1).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',KmidEmidPv};

    pingM.connections(2).direction = [pingM.populations(2).name, '->', pingM.populations(1).name];
    pingM.connections(2).source = pingM.populations(2).name;
    pingM.connections(2).target = pingM.populations(1).name;
    pingM.connections(2).mechanism_list = {'iGABActx'};
    pingM.connections(2).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',KmidPvmidE};

    % Structures: DEEP LAYER (5-6)
    pingD=[];

    % E-cells
    pingD.populations(1).name = ['supE', populationName];
    pingD.populations(1).size = NeSuperficial;
    pingD.populations(1).equations = eqns;
    pingD.populations(1).mechanism_list = cell_type;
    pingD.populations(1).parameters = {'Iapp', 4,'noise', NoiseRate*2};

    % I-cells-SOM
    pingD.populations(2).name = ['supISOM', populationName];
    pingD.populations(2).size = NSomSuperficial;
    pingD.populations(2).equations = eqns;
    pingD.populations(2).mechanism_list = cell_type;
    pingD.populations(2).parameters = {'Iapp',0,'noise', NoiseRate};

    % I-cells-PV
    pingD.populations(3).name = ['supIPV', populationName];
    pingD.populations(3).size = NPvSuperficial;
    pingD.populations(3).equations = eqns;
    pingD.populations(3).mechanism_list = cell_type;
    pingD.populations(3).parameters = {'Iapp',0,'noise', NoiseRate};
    
    % Interlayer connections - Superficial
    pingD.connections(1).direction = [pingD.populations(1).name, '->', pingD.populations(2).name];
    pingD.connections(1).source = pingD.populations(1).name;
    pingD.connections(1).target = pingD.populations(2).name;
    pingD.connections(1).mechanism_list = {'iAMPActx'};
    pingD.connections(1).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',KsupEsupSom};

    pingD.connections(2).direction = [pingD.populations(1).name, '->', pingD.populations(3).name];
    pingD.connections(2).source = pingD.populations(1).name;
    pingD.connections(2).target = pingD.populations(3).name;
    pingD.connections(2).mechanism_list = {'iAMPActx'};
    pingD.connections(2).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',KsupEsupPv};

    pingD.connections(3).direction = [pingD.populations(1).name, '->', pingD.populations(2).name];
    pingD.connections(3).source = pingD.populations(1).name;
    pingD.connections(3).target = pingD.populations(2).name;
    pingD.connections(3).mechanism_list = {'iNMDActx'};
    pingD.connections(3).parameters = {'netcon',KsupEsupSom};

    pingD.connections(4).direction = [pingD.populations(1).name, '->', pingD.populations(3).name];
    pingD.connections(4).source = pingD.populations(1).name;
    pingD.connections(4).target = pingD.populations(3).name;
    pingD.connections(4).mechanism_list = {'iNMDActx'};
    pingD.connections(4).parameters = {'netcon',KsupEsupPv};

    pingD.connections(5).direction = [pingD.populations(2).name, '->', pingD.populations(1).name];
    pingD.connections(5).source = pingD.populations(2).name;
    pingD.connections(5).target = pingD.populations(1).name;
    pingD.connections(5).mechanism_list = {'iGABActx'};
    pingD.connections(5).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',KsupSomsupE};

    pingD.connections(6).direction = [pingD.populations(2).name, '->', pingD.populations(3).name];
    pingD.connections(6).source = pingD.populations(2).name;
    pingD.connections(6).target = pingD.populations(3).name;
    pingD.connections(6).mechanism_list = {'iGABActx'};
    pingD.connections(6).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',KsupSomsupPv};

    pingD.connections(7).direction = [pingD.populations(3).name, '->', pingD.populations(1).name];
    pingD.connections(7).source = pingD.populations(3).name;
    pingD.connections(7).target = pingD.populations(1).name;
    pingD.connections(7).mechanism_list = {'iGABActx'};
    pingD.connections(7).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',KsupPvsupE};

    pingD.connections(8).direction = [pingD.populations(3).name, '->', pingD.populations(2).name];
    pingD.connections(8).source = pingD.populations(3).name;
    pingD.connections(8).target = pingD.populations(2).name;
    pingD.connections(8).mechanism_list = {'iGABActx'};
    pingD.connections(8).parameters = {'gGABAa',gGABAa_ie,'tauGABA',tauGABA_gamma,'netcon',KsupPvsupSom};



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