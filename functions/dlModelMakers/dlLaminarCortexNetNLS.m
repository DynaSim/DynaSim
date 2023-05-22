function y = dlLaminarCortexNetNLS(ModelParameters, populationName)

    NeSuperficial = ModelParameters.NeSuperficial;
    NPvSuperficial = ModelParameters.NPvSuperficial;
    NSomSuperficial = ModelParameters.NSomSuperficial;
    NVipSuperficial = ModelParameters.NVipSuperficial;

    NeMid = ModelParameters.NeMid;
    NPvMid = ModelParameters.NPvMid;
    NSomMid = ModelParameters.NSomMid;
    NVipMid = ModelParameters.NVipMid;

    NeDeep = ModelParameters.NeDeep;
    NPvDeep = ModelParameters.NPvDeep;
    NSomDeep = ModelParameters.NSomDeep;
    NVipDeep = ModelParameters.NVipDeep;

    Nin = ModelParameters.Nin;
    Nout = ModelParameters.Nout;
    Nstim = ModelParameters.Nstim;
    NoiseRate = ModelParameters.NoiseRate;

    totalNeuronCount = NeDeep + NeMid + NeSuperficial + NPvSuperficial + NPvMid + NPvDeep + NSomDeep + NSomMid + NSomSuperficial + NVipSuperficial + NVipMid + NVipDeep;

    fprintf("\n->Initialization of dlLaminarCortexNL Model:\n-->Total neuron count = %d ", totalNeuronCount);
    fprintf("\n-->Superficial (L1-3) Neuron counts: \n--->PYR(E) = %d , PV(I) = %d , SOM(I) = %d , VIP(I) = %d", NeSuperficial, NPvSuperficial, NSomSuperficial, NVipSuperficial);
    fprintf("\n-->Middle (L4) Neuron counts: \n--->PYR(E) = %d , PV(I) = %d , SOM(I) = %d , VIP(I) = %d", NeMid, NPvMid, NSomMid, NVipMid);
    fprintf("\n-->Deep (L5-6) Neuron counts: \n--->PYR(E) = %d , PV(I) = %d , SOM(I) = %d , VIP(I) = %d", NeDeep, NPvDeep, NSomDeep, NVipDeep);

    fprintf("\n--->(E):Excitatory, (I):Inhibitory \n-->Input connections count (terminal) size = %d ", Nin); % Inputs / Stimuli
    fprintf("\n-->Output connections count (terminal) size = %d ", Nout); % Outputs / Probes
    fprintf("\n-->Overall noise rate = %.4f", NoiseRate); % Randomness / Stochasticity
    fprintf("\n->Population name is %s", populationName); % Name tag or suffix for all names of this dsModel

    k1 = 0.15; % Diff. for normal weights (uniform random)
    k2 = 0.45; % Min connectivity weight

    NeAvg = NeSuperficial + NeMid + NeDeep;
    ModelName = populationName;
    populationName = char("_" + populationName);

    % Connectivity matrices

    fprintf("\n***> WARNING! CONNECTIVITY MATRICES NEED TO BE RE-DEFINED.");

    KEESDS = 0.9*ones(NeSuperficial, NeSuperficial);
    KEESDM = 0.9*ones(NeMid, NeMid);
    KEESDD = 0.9*ones(NeDeep, NeDeep);

    % sE->sIsom
    KSupESupSom = k1*rand(NeSuperficial, NSomSuperficial) + 0.7;
%     % sE->sIvip ?
%     KSupESupVip = k1*rand(NeSuperficial, NVipSuperficial) + 0.7;
    % sE->sIpv
    KSupESupPv = k1*rand(NeSuperficial, NPvSuperficial) + 0.7;
    % sE->mE
    KSupEMidE = k1*rand(NeSuperficial, NeMid) + k2;

    % sIsom->sE
    KSupSomSupE = k1*rand(NSomSuperficial, NeSuperficial) + k2;
    % sIsom->sIpv
    KSupSomSupPv = k1*rand(NSomSuperficial, NPvSuperficial) + k2;
    % sIsom->sIvip
    KSupSomSupVip = k1*rand(NSomSuperficial, NVipSuperficial) + k2;

    % sIvip->sIsom
    KSupVipSupSom = k1*rand(NVipSuperficial, NSomSuperficial) + k2;
    % sIvip->sIpv
    KSupVipSupPv = k1*rand(NVipSuperficial, NPvSuperficial) + k2;

    % sIpv->sE
    KSupPvSupE = k1*rand(NPvSuperficial, NeSuperficial) + 0.7;
    % sIpv->sIsom
    KSupPvSupSom = k1*rand(NPvSuperficial, NSomSuperficial) + k2;
    % sIpv->sIsom
    KSupPvSupVip = k1*rand(NPvSuperficial, NVipSuperficial) + k2;

    % mE->mIsom
    KMidEMidSom = k1*rand(NeMid, NSomMid) + 0.7;
%     % mE->mIvip ?
%     KMidEMidVip = k1*rand(NeMid, NVipMid) + 0.7;
    % mE->mIpv
    KMidEMidPv = k1*rand(NeMid, NPvMid) + 0.7;
    % mE->dE
    KMidEDeepE = k1*rand(NeMid, NeDeep) + k2;
    % mE->sE
    KMidESupE = k1*rand(NeMid, NeSuperficial) + k2;
    % dE->sE
    KDeepESupE = k1*rand(NeDeep, NeSuperficial) + k2;

    % mIsom->mE
    KMidSomMidE = k1*rand(NSomMid, NeMid) + k2;
    % mIsom->mIpv
    KMidSomMidPv = k1*rand(NSomMid, NPvMid) + k2;
    % mIsom->mIvip
    KMidSomMidVip = k1*rand(NSomMid, NVipMid) + k2;

    % mIvip->mIsom
    KMidVipMidSom = k1*rand(NVipMid, NSomMid) + k2;
    % mIvip->mIpv
    KMidVipMidPv = k1*rand(NVipMid, NPvMid) + k2;

    % mIpv->mE
    KMidPvMidE = k1*rand(NPvMid, NeMid) + 0.7;
    % mIpv->mIsom
    KMidPvMidSom = k1*rand(NPvMid, NSomMid) + k2;
    % mIpv->mIsom
    KMidPvMidVip = k1*rand(NPvMid, NVipMid) + k2;

    % dE->dIsom
    KDeepEDeepSom = k1*rand(NeDeep, NSomDeep) + 0.7;
%     % dE->dIvip ?
%     KDeepEDeepVip = k1*rand(NeDeep, NVipDeep) + 0.7;
    % dE->dIpv
    KDeepEDeepPv = k1*rand(NeDeep, NPvDeep) + 0.7;
    % dE->mE
    KDeepEMidE = k1*rand(NeDeep, NeMid) + k2;

    % dIsom->dE
    KDeepSomDeepE = k1*rand(NSomDeep, NeDeep) + k2;
    % dIsom->dIpv
    KDeepSomDeepPv = k1*rand(NSomDeep, NPvDeep) + k2;
    % dIsom->dIvip
    KDeepSomDeepVip = k1*rand(NSomDeep, NVipDeep) + k2;

    % dIvip->dIsom
    KDeepVipDeepSom = k1*rand(NVipDeep, NSomDeep) + k2;
    % dIvip->dIpv
    KDeepVipDeepPv = k1*rand(NVipDeep, NPvDeep) + k2;

    % dIpv->dE
    KDeepPvDeepE = k1*rand(NPvDeep, NeDeep) + 0.7;
    % dIpv->dIsom
    KDeepPvDeepSom = k1*rand(NPvDeep, NSomDeep) + k2;
    % dIpv->dIsom
    KDeepPvDeepVip = k1*rand(NPvDeep, NVipDeep) + k2;
    
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
    tauGABA_PvMid = 16.66; % ms, decay time constant of inhibition for gamma (around 50Hz)
    tauGABA_Som = 33.33; % ms, decay time constant of inhibition for beta (around 16Hz)
    tauGABA_PvSup = 16.66;
    tauGABA_PvDeep = 16.66;
    tauGABA_Vip = 166.66;

    tauAMPA_E = 33.33;

    tauAMPA_SM = 33.33; % ms, decay time constant of fast excitation (AMPA)
    tauAMPA_MS = 33.33;
    tauAMPA_MD = 33.33;
    tauAMPA_DM = 33.33;
    
    gBase = 7/(NeAvg^0.5);

    % Maximal synaptic strengths
    gAMPA_MD = .009*gBase; % E->E
    gGABA_SE = .011*gBase; % (SOM)->E
    gGABA_PE = .011*gBase; % (PV)->E
    gGABA_II = .011*gBase; % I->I within layer

    gAMPA_MS = .011*gBase;
    gAMPA_EI_Sup = .004*gBase;
    gAMPA_EI_Mid = .004*gBase;
    gAMPA_EI_Deep = .004*gBase;

    % neuronal dynamics
    eqns = 'dV/dt = (Iapp + @current + noise*randn(1, Npop))/C; Iapp=0; noise=0; C=1; V(0) = -65 + rand(1, Npop)*0;';
    % eqns_apical = 'dV/dt = (Iapp + Icable + @current + noise*randn(1, Npop))/C; Iapp=0; noise=0; C=1; V(0) = -65 + rand(1, Npop)*0;';
    %eqns_stimuli = 'dinp/dt = (rand(1) + 4.5)*(20*(exp(- (t - t1).^2) - exp(- (t - t2).^2)) + noise*randn(1, Npop))/C; f1=4; t1=10; t2=100; noise=0; C=1; inp(0) = -60 - rand(1, Npop)*0;';
    eqns_stimuli = 'dinp/dt = (noise*rand(1) + 1)*(1*(exp(- (t - t1).^2) - exp(- (t - t2).^2)) + noise*randn(1, Npop)); t1=100; t2=200; noise=0; inp(0) = rand(1, Npop)*0;';
    input_amp = 1; % scales the strength of stimulus
    
    % cell type
    ctx_cells = {'iNa','iK', 'iCa', 'iL'};
    cell_type = ctx_cells;

    % BUILDING
    % Structures: SupERFICIAL LAYER (1-2-3)
    SupL=[];

    % E-cells-PYR-Soma
    SupL.populations(1).name = ['SupES', populationName];
    SupL.populations(1).size = NeSuperficial;
    SupL.populations(1).equations = eqns;
    SupL.populations(1).mechanism_list = cell_type;
    SupL.populations(1).parameters = {'Iapp', 0,'noise', NoiseRate*2};

    % E-cells-PYR-Dendrite
    SupL.populations(2).name = ['SupED', populationName];
    SupL.populations(2).size = NeSuperficial;
    SupL.populations(2).equations = eqns;
    SupL.populations(2).mechanism_list = cell_type;
    SupL.populations(2).parameters = {'Iapp', 0,'noise', NoiseRate};

    % I-cells-PV
    SupL.populations(3).name = ['SupIPV', populationName];
    SupL.populations(3).size = NPvSuperficial;
    SupL.populations(3).equations = eqns;
    SupL.populations(3).mechanism_list = cell_type;
    SupL.populations(3).parameters = {'Iapp', 0,'noise', NoiseRate};

    % I-cells-SOM
    SupL.populations(4).name = ['SupISOM', populationName];
    SupL.populations(4).size = NSomSuperficial;
    SupL.populations(4).equations = eqns;
    SupL.populations(4).mechanism_list = cell_type;
    SupL.populations(4).parameters = {'Iapp', 0,'noise', NoiseRate};

    % I-cells-VIP
    SupL.populations(5).name = ['SupIVIP', populationName];
    SupL.populations(5).size = NVipSuperficial;
    SupL.populations(5).equations = eqns;
    SupL.populations(5).mechanism_list = cell_type;
    SupL.populations(5).parameters = {'Iapp', 0,'noise', NoiseRate};
    
    % Interlayer connections - Superficial
    SupL.connections(1).direction = [SupL.populations(1).name, '->', SupL.populations(3).name];
    SupL.connections(1).source = SupL.populations(1).name;
    SupL.connections(1).target = SupL.populations(3).name;
    SupL.connections(1).mechanism_list = {'iAMPActx'};
    SupL.connections(1).parameters = {'gAMPA',gAMPA_EI_Sup,'tauAMPA',tauAMPA_E,'netcon',KSupESupPv};

    SupL.connections(2).direction = [SupL.populations(1).name, '->', SupL.populations(4).name];
    SupL.connections(2).source = SupL.populations(1).name;
    SupL.connections(2).target = SupL.populations(4).name;
    SupL.connections(2).mechanism_list = {'iAMPActx'};
    SupL.connections(2).parameters = {'gAMPA',gAMPA_EI_Sup,'tauAMPA',tauAMPA_E,'netcon',KSupESupSom};

    SupL.connections(3).direction = [SupL.populations(1).name, '->', SupL.populations(3).name];
    SupL.connections(3).source = SupL.populations(1).name;
    SupL.connections(3).target = SupL.populations(3).name;
    SupL.connections(3).mechanism_list = {'iNMDActx'};
    SupL.connections(3).parameters = {'netcon',KSupESupPv};
 
    SupL.connections(4).direction = [SupL.populations(1).name, '->', SupL.populations(4).name];
    SupL.connections(4).source = SupL.populations(1).name;
    SupL.connections(4).target = SupL.populations(4).name;
    SupL.connections(4).mechanism_list = {'iNMDActx'};
    SupL.connections(4).parameters = {'netcon',KSupESupSom};

    SupL.connections(end+1).direction = [SupL.populations(3).name, '->', SupL.populations(1).name];
    SupL.connections(end).source = SupL.populations(2).name;
    SupL.connections(end).target = SupL.populations(1).name;
    SupL.connections(end).mechanism_list = {'iGABActx'};
    SupL.connections(end).parameters = {'gGABAa',gGABA_SE,'tauGABA',tauGABA_PvSup,'netcon',KSupPvSupE};

    SupL.connections(end+1).direction = [SupL.populations(3).name, '->', SupL.populations(4).name];
    SupL.connections(end).source = SupL.populations(3).name;
    SupL.connections(end).target = SupL.populations(4).name;
    SupL.connections(end).mechanism_list = {'iGABActx'};
    SupL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_PvDeep,'netcon',KSupPvSupSom};

    SupL.connections(end+1).direction = [SupL.populations(4).name, '->', SupL.populations(2).name];
    SupL.connections(end).source = SupL.populations(4).name;
    SupL.connections(end).target = SupL.populations(2).name;
    SupL.connections(end).mechanism_list = {'iGABActx'};
    SupL.connections(end).parameters = {'gGABAa',gGABA_PE,'tauGABA',tauGABA_Som,'netcon',KSupSomSupE};

    SupL.connections(end+1).direction = [SupL.populations(4).name, '->', SupL.populations(3).name];
    SupL.connections(end).source = SupL.populations(4).name;
    SupL.connections(end).target = SupL.populations(3).name;
    SupL.connections(end).mechanism_list = {'iGABActx'};
    SupL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_Som,'netcon',KSupSomSupPv};

    SupL.connections(end+1).direction = [SupL.populations(5).name, '->', SupL.populations(3).name];
    SupL.connections(end).source = SupL.populations(5).name;
    SupL.connections(end).target = SupL.populations(3).name;
    SupL.connections(end).mechanism_list = {'iGABActx'};
    SupL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_Vip,'netcon',KSupVipSupPv};

    SupL.connections(end+1).direction = [SupL.populations(5).name, '->', SupL.populations(4).name];
    SupL.connections(end).source = SupL.populations(5).name;
    SupL.connections(end).target = SupL.populations(4).name;
    SupL.connections(end).mechanism_list = {'iGABActx'};
    SupL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_Vip,'netcon',KSupVipSupSom};

    SupL.connections(end+1).direction = [SupL.populations(3).name, '->', SupL.populations(5).name];
    SupL.connections(end).source = SupL.populations(3).name;
    SupL.connections(end).target = SupL.populations(5).name;
    SupL.connections(end).mechanism_list = {'iGABActx'};
    SupL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_PvDeep,'netcon',KSupPvSupVip};

    SupL.connections(end+1).direction = [SupL.populations(4).name, '->', SupL.populations(5).name];
    SupL.connections(end).source = SupL.populations(4).name;
    SupL.connections(end).target = SupL.populations(5).name;
    SupL.connections(end).mechanism_list = {'iGABActx'};
    SupL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_Som,'netcon',KSupSomSupVip};

    SupL.connections(end+1).direction = [SupL.populations(1).name, '->', SupL.populations(2).name];
    SupL.connections(end).source = SupL.populations(1).name;
    SupL.connections(end).target = SupL.populations(2).name;
    SupL.connections(end).mechanism_list = {'iAMPActx'};
    SupL.connections(end).parameters = {'gGABAa',gAMPA_MS,'tauGABA',tauAMPA_SM,'netcon',KEESDS};

    % Structures: MIDDLE LAYER (4)
    MidL=[];

    % E-cells-PYR-Soma
    MidL.populations(1).name = ['MidES', populationName];
    MidL.populations(1).size = NeMid;
    MidL.populations(1).equations = eqns;
    MidL.populations(1).mechanism_list = cell_type;
    MidL.populations(1).parameters = {'Iapp', 0,'noise', NoiseRate*2};

    % E-cells-PYR-Dendrite
    MidL.populations(2).name = ['MidED', populationName];
    MidL.populations(2).size = NeMid;
    MidL.populations(2).equations = eqns;
    MidL.populations(2).mechanism_list = cell_type;
    MidL.populations(2).parameters = {'Iapp', 0,'noise', NoiseRate};

    % I-cells-PV
    MidL.populations(3).name = ['MidIPV', populationName];
    MidL.populations(3).size = NPvMid;
    MidL.populations(3).equations = eqns;
    MidL.populations(3).mechanism_list = cell_type;
    MidL.populations(3).parameters = {'Iapp', 0,'noise', NoiseRate};

    % I-cells-SOM
    MidL.populations(4).name = ['MidISOM', populationName];
    MidL.populations(4).size = NSomMid;
    MidL.populations(4).equations = eqns;
    MidL.populations(4).mechanism_list = cell_type;
    MidL.populations(4).parameters = {'Iapp', 0,'noise', NoiseRate};

    % I-cells-VIP
    MidL.populations(5).name = ['MidIVIP', populationName];
    MidL.populations(5).size = NVipMid;
    MidL.populations(5).equations = eqns;
    MidL.populations(5).mechanism_list = cell_type;
    MidL.populations(5).parameters = {'Iapp', 0,'noise', NoiseRate};
    
    % Interlayer connections - Middle
    MidL.connections(1).direction = [MidL.populations(1).name, '->', MidL.populations(3).name];
    MidL.connections(1).source = MidL.populations(1).name;
    MidL.connections(1).target = MidL.populations(3).name;
    MidL.connections(1).mechanism_list = {'iAMPActx'};
    MidL.connections(1).parameters = {'gAMPA',gAMPA_EI_Mid,'tauAMPA',tauAMPA_E,'netcon',KMidEMidPv};

    MidL.connections(2).direction = [MidL.populations(1).name, '->', MidL.populations(4).name];
    MidL.connections(2).source = MidL.populations(1).name;
    MidL.connections(2).target = MidL.populations(4).name;
    MidL.connections(2).mechanism_list = {'iAMPActx'};
    MidL.connections(2).parameters = {'gAMPA',gAMPA_EI_Mid,'tauAMPA',tauAMPA_E,'netcon',KMidEMidSom};

    MidL.connections(3).direction = [MidL.populations(1).name, '->', MidL.populations(3).name];
    MidL.connections(3).source = MidL.populations(1).name;
    MidL.connections(3).target = MidL.populations(2).name;
    MidL.connections(3).mechanism_list = {'iNMDActx'};
    MidL.connections(3).parameters = {'netcon',KMidEMidPv};
 
    MidL.connections(4).direction = [MidL.populations(1).name, '->', MidL.populations(4).name];
    MidL.connections(4).source = MidL.populations(1).name;
    MidL.connections(4).target = MidL.populations(3).name;
    MidL.connections(4).mechanism_list = {'iNMDActx'};
    MidL.connections(4).parameters = {'netcon',KMidEMidSom};

    MidL.connections(end+1).direction = [MidL.populations(3).name, '->', MidL.populations(1).name];
    MidL.connections(end).source = MidL.populations(2).name;
    MidL.connections(end).target = MidL.populations(1).name;
    MidL.connections(end).mechanism_list = {'iGABActx'};
    MidL.connections(end).parameters = {'gGABAa',gGABA_SE,'tauGABA',tauGABA_PvMid,'netcon',KMidPvMidE};

    MidL.connections(end+1).direction = [MidL.populations(3).name, '->', MidL.populations(4).name];
    MidL.connections(end).source = MidL.populations(3).name;
    MidL.connections(end).target = MidL.populations(4).name;
    MidL.connections(end).mechanism_list = {'iGABActx'};
    MidL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_PvDeep,'netcon',KMidPvMidSom};

    MidL.connections(end+1).direction = [MidL.populations(4).name, '->', MidL.populations(2).name];
    MidL.connections(end).source = MidL.populations(4).name;
    MidL.connections(end).target = MidL.populations(2).name;
    MidL.connections(end).mechanism_list = {'iGABActx'};
    MidL.connections(end).parameters = {'gGABAa',gGABA_PE,'tauGABA',tauGABA_Som,'netcon',KMidSomMidE};

    MidL.connections(end+1).direction = [MidL.populations(4).name, '->', MidL.populations(3).name];
    MidL.connections(end).source = MidL.populations(4).name;
    MidL.connections(end).target = MidL.populations(3).name;
    MidL.connections(end).mechanism_list = {'iGABActx'};
    MidL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_Som,'netcon',KMidSomMidPv};

    MidL.connections(end+1).direction = [MidL.populations(5).name, '->', MidL.populations(3).name];
    MidL.connections(end).source = MidL.populations(5).name;
    MidL.connections(end).target = MidL.populations(3).name;
    MidL.connections(end).mechanism_list = {'iGABActx'};
    MidL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_Vip,'netcon',KMidVipMidPv};

    MidL.connections(end+1).direction = [MidL.populations(5).name, '->', MidL.populations(4).name];
    MidL.connections(end).source = MidL.populations(5).name;
    MidL.connections(end).target = MidL.populations(4).name;
    MidL.connections(end).mechanism_list = {'iGABActx'};
    MidL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_Vip,'netcon',KMidVipMidSom};

    MidL.connections(end+1).direction = [MidL.populations(3).name, '->', MidL.populations(5).name];
    MidL.connections(end).source = MidL.populations(3).name;
    MidL.connections(end).target = MidL.populations(5).name;
    MidL.connections(end).mechanism_list = {'iGABActx'};
    MidL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_PvDeep,'netcon',KMidPvMidVip};

    MidL.connections(end+1).direction = [MidL.populations(4).name, '->', MidL.populations(5).name];
    MidL.connections(end).source = MidL.populations(4).name;
    MidL.connections(end).target = MidL.populations(5).name;
    MidL.connections(end).mechanism_list = {'iGABActx'};
    MidL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_Som,'netcon',KMidSomMidVip};

    MidL.connections(end+1).direction = [MidL.populations(1).name, '->', MidL.populations(2).name];
    MidL.connections(end).source = MidL.populations(1).name;
    MidL.connections(end).target = MidL.populations(2).name;
    MidL.connections(end).mechanism_list = {'iAMPActx'};
    MidL.connections(end).parameters = {'gGABAa',gAMPA_MS,'tauGABA',tauAMPA_SM,'netcon',KEESDM};

    % Structures: DEEP LAYER (5-6)
    DeepL=[];

    % E-cells-PYR-Soma
    DeepL.populations(1).name = ['DeepES', populationName];
    DeepL.populations(1).size = NeDeep;
    DeepL.populations(1).equations = eqns;
    DeepL.populations(1).mechanism_list = cell_type;
    DeepL.populations(1).parameters = {'Iapp', 0,'noise', NoiseRate*2};

    % E-cells-PYR-Dendrite
    DeepL.populations(2).name = ['DeepED', populationName];
    DeepL.populations(2).size = NeDeep;
    DeepL.populations(2).equations = eqns;
    DeepL.populations(2).mechanism_list = cell_type;
    DeepL.populations(2).parameters = {'Iapp', 0,'noise', NoiseRate};

    % I-cells-PV
    DeepL.populations(3).name = ['DeepIPV', populationName];
    DeepL.populations(3).size = NPvDeep;
    DeepL.populations(3).equations = eqns;
    DeepL.populations(3).mechanism_list = cell_type;
    DeepL.populations(3).parameters = {'Iapp', 0,'noise', NoiseRate};

    % I-cells-SOM
    DeepL.populations(4).name = ['DeepISOM', populationName];
    DeepL.populations(4).size = NSomDeep;
    DeepL.populations(4).equations = eqns;
    DeepL.populations(4).mechanism_list = cell_type;
    DeepL.populations(4).parameters = {'Iapp', 0,'noise', NoiseRate};

    % I-cells-VIP
    DeepL.populations(5).name = ['DeepIVIP', populationName];
    DeepL.populations(5).size = NVipDeep;
    DeepL.populations(5).equations = eqns;
    DeepL.populations(5).mechanism_list = cell_type;
    DeepL.populations(5).parameters = {'Iapp', 0,'noise', NoiseRate};
    
    % Interlayer connections - Deepdle
    DeepL.connections(1).direction = [DeepL.populations(1).name, '->', DeepL.populations(3).name];
    DeepL.connections(1).source = DeepL.populations(1).name;
    DeepL.connections(1).target = DeepL.populations(3).name;
    DeepL.connections(1).mechanism_list = {'iAMPActx'};
    DeepL.connections(1).parameters = {'gAMPA',gAMPA_EI_Deep,'tauAMPA',tauAMPA_E,'netcon',KDeepEDeepPv};

    DeepL.connections(2).direction = [DeepL.populations(1).name, '->', DeepL.populations(4).name];
    DeepL.connections(2).source = DeepL.populations(1).name;
    DeepL.connections(2).target = DeepL.populations(4).name;
    DeepL.connections(2).mechanism_list = {'iAMPActx'};
    DeepL.connections(2).parameters = {'gAMPA',gAMPA_EI_Deep,'tauAMPA',tauAMPA_E,'netcon',KDeepEDeepSom};

    DeepL.connections(3).direction = [DeepL.populations(1).name, '->', DeepL.populations(3).name];
    DeepL.connections(3).source = DeepL.populations(1).name;
    DeepL.connections(3).target = DeepL.populations(2).name;
    DeepL.connections(3).mechanism_list = {'iNMDActx'};
    DeepL.connections(3).parameters = {'netcon',KDeepEDeepPv};
 
    DeepL.connections(4).direction = [DeepL.populations(1).name, '->', DeepL.populations(4).name];
    DeepL.connections(4).source = DeepL.populations(1).name;
    DeepL.connections(4).target = DeepL.populations(3).name;
    DeepL.connections(4).mechanism_list = {'iNMDActx'};
    DeepL.connections(4).parameters = {'netcon',KDeepEDeepSom};

    DeepL.connections(end+1).direction = [DeepL.populations(3).name, '->', DeepL.populations(1).name];
    DeepL.connections(end).source = DeepL.populations(2).name;
    DeepL.connections(end).target = DeepL.populations(1).name;
    DeepL.connections(end).mechanism_list = {'iGABActx'};
    DeepL.connections(end).parameters = {'gGABAa',gGABA_SE,'tauGABA',tauGABA_PvDeep,'netcon',KDeepPvDeepE};

    DeepL.connections(end+1).direction = [DeepL.populations(3).name, '->', DeepL.populations(4).name];
    DeepL.connections(end).source = DeepL.populations(3).name;
    DeepL.connections(end).target = DeepL.populations(4).name;
    DeepL.connections(end).mechanism_list = {'iGABActx'};
    DeepL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_PvDeep,'netcon',KDeepPvDeepSom};

    DeepL.connections(end+1).direction = [DeepL.populations(4).name, '->', DeepL.populations(2).name];
    DeepL.connections(end).source = DeepL.populations(4).name;
    DeepL.connections(end).target = DeepL.populations(2).name;
    DeepL.connections(end).mechanism_list = {'iGABActx'};
    DeepL.connections(end).parameters = {'gGABAa',gGABA_PE,'tauGABA',tauGABA_Som,'netcon',KDeepSomDeepE};

    DeepL.connections(end+1).direction = [DeepL.populations(4).name, '->', DeepL.populations(3).name];
    DeepL.connections(end).source = DeepL.populations(4).name;
    DeepL.connections(end).target = DeepL.populations(3).name;
    DeepL.connections(end).mechanism_list = {'iGABActx'};
    DeepL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_Som,'netcon',KDeepSomDeepPv};

    DeepL.connections(end+1).direction = [DeepL.populations(5).name, '->', DeepL.populations(3).name];
    DeepL.connections(end).source = DeepL.populations(5).name;
    DeepL.connections(end).target = DeepL.populations(3).name;
    DeepL.connections(end).mechanism_list = {'iGABActx'};
    DeepL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_Vip,'netcon',KDeepVipDeepPv};

    DeepL.connections(end+1).direction = [DeepL.populations(5).name, '->', DeepL.populations(4).name];
    DeepL.connections(end).source = DeepL.populations(5).name;
    DeepL.connections(end).target = DeepL.populations(4).name;
    DeepL.connections(end).mechanism_list = {'iGABActx'};
    DeepL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_Vip,'netcon',KDeepVipDeepSom};

    DeepL.connections(end+1).direction = [DeepL.populations(3).name, '->', DeepL.populations(5).name];
    DeepL.connections(end).source = DeepL.populations(3).name;
    DeepL.connections(end).target = DeepL.populations(5).name;
    DeepL.connections(end).mechanism_list = {'iGABActx'};
    DeepL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_PvDeep,'netcon',KDeepPvDeepVip};

    DeepL.connections(end+1).direction = [DeepL.populations(4).name, '->', DeepL.populations(5).name];
    DeepL.connections(end).source = DeepL.populations(4).name;
    DeepL.connections(end).target = DeepL.populations(5).name;
    DeepL.connections(end).mechanism_list = {'iGABActx'};
    DeepL.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_Som,'netcon',KDeepSomDeepVip};

    DeepL.connections(end+1).direction = [DeepL.populations(1).name, '->', DeepL.populations(2).name];
    DeepL.connections(end).source = DeepL.populations(1).name;
    DeepL.connections(end).target = DeepL.populations(2).name;
    DeepL.connections(end).mechanism_list = {'iAMPActx'};
    DeepL.connections(end).parameters = {'gGABAa',gAMPA_MS,'tauGABA',tauAMPA_SM,'netcon',KEESDD};

    % PING template
    IOping = cell(Nstim, 1);

    for i = 1:Nstim

        % E-cells
        IOping{i}.populations(1).name = ['PreStimuli', char(64+i), populationName];
        IOping{i}.populations(1).size = Nin;
        IOping{i}.populations(1).equations = eqns_stimuli;
        IOping{i}.populations(1).mechanism_list = cell_type;
        IOping{i}.populations(1).parameters = {'f1', 1,'noise', 0};
    
        % I-cells
        IOping{i}.populations(2).name = ['PostStimuli', char(64+i), populationName];
        IOping{i}.populations(2).size = Nin;
        IOping{i}.populations(2).equations = eqns_stimuli;
        IOping{i}.populations(2).mechanism_list = cell_type;
        IOping{i}.populations(2).parameters = {'f1', 1,'noise', 0};

        IOping{i}.connections(1).direction = [IOping{i}.populations(1).name, '->', IOping{i}.populations(2).name];
        IOping{i}.connections(1).source = IOping{i}.populations(1).name;
        IOping{i}.connections(1).target = IOping{i}.populations(2).name;
        IOping{i}.connections(1).mechanism_list = {'iAMPActx'};
        IOping{i}.connections(1).parameters = {'gAMPA', 1e-7,'tauGABA', 1e+4,'netcon', KnullIO};
        
    end

    s = dsCombineSpecifications(SupL, MidL, DeepL, IOping{1}, IOping{2}, IOping{3});

    % connect the layers and inputs
    fprintf("\n-->Connecting separate layers and inputs:");
    tempconn = zeros(Nin, NeMid);

    % Inputs/stimuli -> midE

    for k = 1:Nstim

        Aconn = tempconn;
        Aconn(:, subLayerIndicesInM(1, k):subLayerIndicesInM(2, k)) =  0.74;
    
        c = length(s.connections) + 1;
        s.connections(c).direction = [IOping{k}.populations(1).name, '->', MidL.populations(1).name];
        s.connections(c).source = IOping{k}.populations(1).name;
        s.connections(c).target = MidL.populations(1).name;
        s.connections(c).mechanism_list={'input_connection'};
        s.connections(c).parameters={'input_amp',input_amp,'netcon',Aconn};
    
        c = length(s.connections) + 1;
        s.connections(c).direction = [IOping{k}.populations(2).name, '->', MidL.populations(1).name];
        s.connections(c).source = IOping{k}.populations(2).name;
        s.connections(c).target = MidL.populations(1).name;
        s.connections(c).mechanism_list={'input_connection'};
        s.connections(c).parameters={'input_amp',input_amp,'netcon',Aconn};

    end
    
    % SupEmidE
    c = length(s.connections)+1;
    s.connections(c).direction = [SupL.populations(2).name, '->', MidL.populations(1).name];
    s.connections(c).source = SupL.populations(2).name;
    s.connections(c).target = MidL.populations(1).name;
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_MS,'tauAMPA',tauAMPA_SM,'netcon',KSupEMidE};

    % midESupE
    c = length(s.connections)+1;
    s.connections(c).direction = [MidL.populations(2).name, '->', SupL.populations(1).name];
    s.connections(c).source = MidL.populations(2).name;
    s.connections(c).target = SupL.populations(1).name;
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_MS,'tauAMPA',tauAMPA_MS,'netcon',KMidESupE};

    % midEdeepE
    c = length(s.connections)+1;
    s.connections(c).direction = [MidL.populations(2).name, '->', DeepL.populations(1).name];
    s.connections(c).source = MidL.populations(2).name;
    s.connections(c).target = DeepL.populations(1).name;
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_MD,'tauAMPA',tauAMPA_MD,'netcon',KMidEDeepE};
    
    % deepEmidE
    c = length(s.connections)+1;
    s.connections(c).direction = [DeepL.populations(2).name, '->', MidL.populations(1).name];
    s.connections(c).source = DeepL.populations(2).name;
    s.connections(c).target = MidL.populations(1).name;
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_MD,'tauAMPA',tauAMPA_DM,'netcon',KDeepEMidE};

    % deepEsupE
    c = length(s.connections)+1;
    s.connections(c).direction = [DeepL.populations(2).name, '->', SupL.populations(1).name];
    s.connections(c).source = DeepL.populations(2).name;
    s.connections(c).target = SupL.populations(1).name;
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_MD,'tauAMPA',tauAMPA_DM,'netcon',KDeepESupE};

    % Outputs: deepE is the recommended output layer.
    y = s;
    fprintf("\n->Initialization of dsModel ""%s"" is done. \n", ModelName);

end