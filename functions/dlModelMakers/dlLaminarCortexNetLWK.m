function y = dlLaminarCortexNetLWK(ModelParameters, populationName)

    NeSuperficial = ModelParameters.NeSuperficial;
    NSomSuperficial = ModelParameters.NSomSuperficial;
    NPvSuperficial = ModelParameters.NPvSuperficial;

    NeMid = ModelParameters.NeMid;
    NSomMid = 0*ModelParameters.NSomMid;
    NPvMid = ModelParameters.NPvMid;

    NeDeep = ModelParameters.NeDeep;
    NSomDeep = ModelParameters.NSomDeep;
    NPvDeep = ModelParameters.NPvDeep;

    Nin = ModelParameters.Nin;
    Nout = ModelParameters.Nout;

    Nstim = ModelParameters.Nstim;
    NoiseRate = ModelParameters.NoiseRate;

    fprintf("\n>Initialization of dlLaminarCortex Model: ");
    fprintf("\n-->Superficial (L1-3) excitatory neurons count = %d , SOM inhibitory = %d , PV inhibitory = %d ", NeSuperficial, NSomSuperficial, NPvSuperficial); % S
    fprintf("\n-->Middle (L4) excitatory neurons count = %d , SOM inhibitory = %d , PV inhibitory = %d ", NeMid, NSomMid, NPvMid); % M
    fprintf("\n-->Deep (L5-6) excitatory neurons count = %d , SOM inhibitory = %d , PV inhibitory = %d ", NeDeep, NSomDeep, NPvDeep); % D 

    fprintf("\n-->Input connections count (terminal) size = %d ", Nin); % Inputs / Stimuli
    fprintf("\n-->Output connections count (terminal) size = %d ", Nout); % Outputs / Probes
    fprintf("\n-->Overall noise rate = %.4f", NoiseRate); % Randomness / Stochasticity
    fprintf("\n--->Population name is %s", populationName); % Name tag or suffix for all names of this dsModel

    k1 = 0.15; % Diff. for normal weights (uniform random)
    k2 = 0.45; % Min connectivity weight

    NeAvg = NeSuperficial + NeMid + NeDeep;

    populationName = ['x', populationName];
    % Connectivity matrices

    % sE->sIsom
    KsupEsupSom = k1*rand(NeSuperficial, NSomSuperficial) + 0.7;
    % sE->sIpv
    KsupEsupPv = k1*rand(NeSuperficial, NPvSuperficial) + 0.7;
    % sE->dE
    KsupEmidE = k1*rand(NeSuperficial, NeMid) + k2;

    % sIsom->sE
    KsupSomsupE = k1*rand(NSomSuperficial, NeSuperficial) + k2;
    % sIsom->sIpv
    KsupSomsupPv = k1*rand(NSomSuperficial, NPvSuperficial) + k2;

    % sIpv->sE
    KsupPvsupE = k1*rand(NPvSuperficial, NeSuperficial) + 0.7;
    % sIpv->sIsom
    KsupPvsupSom = k1*rand(NPvSuperficial, NSomSuperficial) + k2;

    % mE->mIpv
    KmidEmidPv = k1*rand(NeMid, NPvMid) + k2;
    % mE->sE
    KmidEsupE = k1*rand(NeMid, NeSuperficial) + k2;
    % mE->dE
    KmidEdeepE = k1*rand(NeMid, NeDeep) + k2;

    % mIpv->mE
    KmidPvmidE = k1*rand(NPvMid, NeMid) + k2;

    % dE->dIsom
    KdeepEdeepSom = k1*rand(NeDeep, NSomDeep) + k2;
    % dE->dIpv
    KdeepEdeepPv = k1*rand(NeDeep, NPvDeep) + k2;
    % dE->sIpv
    KdeepEmidE = k1*rand(NeDeep, NeMid) + k2;

    % dIpv->dE
    KdeepPvdeepE = k1*rand(NPvDeep, NeDeep) + k2;
    % dIpv->dISom
    KdeepPvdeepSom = k1*rand(NPvDeep, NSomDeep) + k2;

    % dIsom->dE
    KdeepSomdeepE = k1*rand(NSomDeep, NeDeep) + k2;
    % dIsom->dIpv
    KdeepSomdeepPv = k1*rand(NSomDeep, NPvDeep) + k2;
    
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
    tauGABA_PvMid = 53; % ms, decay time constant of inhibition for gamma (around 50Hz)
    tauGABA_Som = 151; % ms, decay time constant of inhibition for beta (around 25Hz)
    tauGABA_PvSup = 57;
    tauGABA_PvDeep = 91;

    tauAMPA_E = 171;

    tauAMPA_SM = 49; % ms, decay time constant of fast excitation (AMPA)
    tauAMPA_MS = 47;
    tauAMPA_MD = 111;
    tauAMPA_DM = 91;
    
    gBase = 7/(NeAvg^0.5);

    % Maximal synaptic strengths
    gAMPA_MD = .009*gBase; % E->E
    gGABA_SE = .017*gBase; % (SOM)->E
    gGABA_PE = .011*gBase; % (PV)->E
    gGABA_II = .011*gBase; % I->I within layer

    gAMPA_MS = .017*gBase;
    gAMPA_EI_sup = .011*gBase;
    gAMPA_EI_mid = .014*gBase;
    gAMPA_EI_deep = .010*gBase;

    % neuronal dynamics
    eqns = 'dV/dt = (Iapp + @current + noise*randn(1, Npop))/C; Iapp=0; noise=0; C=1; V(0) = -65 + rand(1, Npop)*0;';
    %eqns_stimuli = 'dinp/dt = (rand(1) + 4.5)*(20*(exp(- (t - t1).^2) - exp(- (t - t2).^2)) + noise*randn(1, Npop))/C; f1=4; t1=10; t2=100; noise=0; C=1; inp(0) = -60 - rand(1, Npop)*0;';
    eqns_stimuli = 'dinp/dt = (noise*rand(1) + 1)*(1*(exp(- (t - t1).^2) - exp(- (t - t2).^2)) + noise*randn(1, Npop)); t1=100; t2=200; noise=0; inp(0) = rand(1, Npop)*0;';
    input_amp = 1; % scales the strength of stimulus
    
    % cell type
    ctx_cells = {'iNa','iK'};
    cell_type = ctx_cells;

    % BUILDING
    % Structures: SUPERFICIAL LAYER (1-2-3)
    pingS=[];

    % E-cells-PYR
    pingS.populations(1).name = ['supE', populationName];
    pingS.populations(1).size = NeSuperficial;
    pingS.populations(1).equations = eqns;
    pingS.populations(1).mechanism_list = cell_type;
    pingS.populations(1).parameters = {'Iapp', 0,'noise', NoiseRate*2};

    % I-cells-SOM
    pingS.populations(2).name = ['supISOM', populationName];
    pingS.populations(2).size = NSomSuperficial;
    pingS.populations(2).equations = eqns;
    pingS.populations(2).mechanism_list = cell_type;
    pingS.populations(2).parameters = {'Iapp', 0,'noise', NoiseRate};

    % I-cells-PV
    pingS.populations(3).name = ['supIPV', populationName];
    pingS.populations(3).size = NPvSuperficial;
    pingS.populations(3).equations = eqns;
    pingS.populations(3).mechanism_list = cell_type;
    pingS.populations(3).parameters = {'Iapp', 0,'noise', NoiseRate};
    
    % Interlayer connections - Superficial
    pingS.connections(1).direction = [pingS.populations(1).name, '->', pingS.populations(2).name];
    pingS.connections(1).source = pingS.populations(1).name;
    pingS.connections(1).target = pingS.populations(2).name;
    pingS.connections(1).mechanism_list = {'iAMPActx'};
    pingS.connections(1).parameters = {'gAMPA',gAMPA_EI_sup,'tauAMPA',tauAMPA_E,'netcon',KsupEsupSom};

    pingS.connections(2).direction = [pingS.populations(1).name, '->', pingS.populations(3).name];
    pingS.connections(2).source = pingS.populations(1).name;
    pingS.connections(2).target = pingS.populations(3).name;
    pingS.connections(2).mechanism_list = {'iAMPActx'};
    pingS.connections(2).parameters = {'gAMPA',gAMPA_EI_sup,'tauAMPA',tauAMPA_E,'netcon',KsupEsupPv};

%     pingS.connections(3).direction = [pingS.populations(1).name, '->', pingS.populations(2).name];
%     pingS.connections(3).source = pingS.populations(1).name;
%     pingS.connections(3).target = pingS.populations(2).name;
%     pingS.connections(3).mechanism_list = {'iNMDActx'};
%     pingS.connections(3).parameters = {'netcon',KsupEsupSom};
% 
%     pingS.connections(4).direction = [pingS.populations(1).name, '->', pingS.populations(3).name];
%     pingS.connections(4).source = pingS.populations(1).name;
%     pingS.connections(4).target = pingS.populations(3).name;
%     pingS.connections(4).mechanism_list = {'iNMDActx'};
%     pingS.connections(4).parameters = {'netcon',KsupEsupPv};

    pingS.connections(end+1).direction = [pingS.populations(2).name, '->', pingS.populations(1).name];
    pingS.connections(end).source = pingS.populations(2).name;
    pingS.connections(end).target = pingS.populations(1).name;
    pingS.connections(end).mechanism_list = {'iGABActx'};
    pingS.connections(end).parameters = {'gGABAa',gGABA_SE,'tauGABA',tauGABA_Som,'netcon',KsupSomsupE};

    pingS.connections(end+1).direction = [pingS.populations(2).name, '->', pingS.populations(3).name];
    pingS.connections(end).source = pingS.populations(2).name;
    pingS.connections(end).target = pingS.populations(3).name;
    pingS.connections(end).mechanism_list = {'iGABActx'};
    pingS.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_Som,'netcon',KsupSomsupPv};

    pingS.connections(end+1).direction = [pingS.populations(3).name, '->', pingS.populations(1).name];
    pingS.connections(end).source = pingS.populations(3).name;
    pingS.connections(end).target = pingS.populations(1).name;
    pingS.connections(end).mechanism_list = {'iGABActx'};
    pingS.connections(end).parameters = {'gGABAa',gGABA_PE,'tauGABA',tauGABA_PvSup,'netcon',KsupPvsupE};

    pingS.connections(end+1).direction = [pingS.populations(3).name, '->', pingS.populations(2).name];
    pingS.connections(end).source = pingS.populations(3).name;
    pingS.connections(end).target = pingS.populations(2).name;
    pingS.connections(end).mechanism_list = {'iGABActx'};
    pingS.connections(end).parameters = {'gGABAa',gGABA_II,'tauGABA',tauGABA_PvSup,'netcon',KsupPvsupSom};

    % Structures: MIDDLE LAYER (4)
    pingM=[];

    % E-cells
    pingM.populations(1).name = ['midE', populationName];
    pingM.populations(1).size = NeMid;
    pingM.populations(1).equations = eqns;
    pingM.populations(1).mechanism_list = cell_type;
    pingM.populations(1).parameters = {'Iapp', 0,'noise', NoiseRate*2};

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
    pingM.connections(1).parameters = {'gAMPA',gAMPA_EI_mid,'tauAMPA',tauAMPA_E,'netcon',KmidEmidPv};

    pingM.connections(2).direction = [pingM.populations(2).name, '->', pingM.populations(1).name];
    pingM.connections(2).source = pingM.populations(2).name;
    pingM.connections(2).target = pingM.populations(1).name;
    pingM.connections(2).mechanism_list = {'iGABActx'};
    pingM.connections(2).parameters = {'gGABAa',gAMPA_EI_mid,'tauGABA',tauGABA_PvMid,'netcon',KmidPvmidE};

%     pingIM = pingM;
    % Structures: DEEP LAYER (5-6)
    pingD=[];

    % E-cells
    pingD.populations(1).name = ['deepE', populationName];
    pingD.populations(1).size = NeDeep;
    pingD.populations(1).equations = eqns;
    pingD.populations(1).mechanism_list = cell_type;
    pingD.populations(1).parameters = {'Iapp', 0,'noise', NoiseRate*2};

    % I-cells-SOM
    pingD.populations(2).name = ['deepISOM', populationName];
    pingD.populations(2).size = NSomDeep;
    pingD.populations(2).equations = eqns;
    pingD.populations(2).mechanism_list = cell_type;
    pingD.populations(2).parameters = {'Iapp',0,'noise', NoiseRate};

    % I-cells-PV
    pingD.populations(3).name = ['deepIPV', populationName];
    pingD.populations(3).size = NPvDeep;
    pingD.populations(3).equations = eqns;
    pingD.populations(3).mechanism_list = cell_type;
    pingD.populations(3).parameters = {'Iapp',0,'noise', NoiseRate};
    
    % Interlayer connections - Deep
    pingD.connections(1).direction = [pingD.populations(1).name, '->', pingD.populations(2).name];
    pingD.connections(1).source = pingD.populations(1).name;
    pingD.connections(1).target = pingD.populations(2).name;
    pingD.connections(1).mechanism_list = {'iAMPActx'};
    pingD.connections(1).parameters = {'gAMPA',gAMPA_EI_deep,'tauAMPA',tauAMPA_E*2,'netcon',KdeepEdeepSom};

    pingD.connections(2).direction = [pingD.populations(1).name, '->', pingD.populations(3).name];
    pingD.connections(2).source = pingD.populations(1).name;
    pingD.connections(2).target = pingD.populations(3).name;
    pingD.connections(2).mechanism_list = {'iAMPActx'};
    pingD.connections(2).parameters = {'gAMPA',gAMPA_EI_deep,'tauAMPA',tauAMPA_E*2,'netcon',KdeepEdeepPv};

    pingD.connections(3).direction = [pingD.populations(2).name, '->', pingD.populations(1).name];
    pingD.connections(3).source = pingD.populations(2).name;
    pingD.connections(3).target = pingD.populations(1).name;
    pingD.connections(3).mechanism_list = {'iGABActx'};
    pingD.connections(3).parameters = {'gGABAa',gAMPA_EI_deep,'tauGABA',tauGABA_Som*2,'netcon',KdeepSomdeepE};

    pingD.connections(4).direction = [pingD.populations(2).name, '->', pingD.populations(3).name];
    pingD.connections(4).source = pingD.populations(2).name;
    pingD.connections(4).target = pingD.populations(3).name;
    pingD.connections(4).mechanism_list = {'iGABActx'};
    pingD.connections(4).parameters = {'gGABAa',gAMPA_EI_deep,'tauGABA',tauGABA_Som*2,'netcon',KdeepSomdeepPv};

    pingD.connections(5).direction = [pingD.populations(3).name, '->', pingD.populations(1).name];
    pingD.connections(5).source = pingD.populations(3).name;
    pingD.connections(5).target = pingD.populations(1).name;
    pingD.connections(5).mechanism_list = {'iGABActx'};
    pingD.connections(5).parameters = {'gGABAa',gAMPA_EI_deep,'tauGABA',tauGABA_PvDeep,'netcon',KdeepPvdeepE};

    pingD.connections(6).direction = [pingD.populations(3).name, '->', pingD.populations(2).name];
    pingD.connections(6).source = pingD.populations(3).name;
    pingD.connections(6).target = pingD.populations(2).name;
    pingD.connections(6).mechanism_list = {'iGABActx'};
    pingD.connections(6).parameters = {'gGABAa',gAMPA_EI_deep,'tauGABA',tauGABA_PvDeep,'netcon',KdeepPvdeepSom};

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

    s = dsCombineSpecifications(pingS, pingM, pingD, IOping{1}, IOping{2}, IOping{3});

    % connect the layers and inputs
    fprintf("\n--->Connecting separate layers and inputs:");
    tempconn = zeros(Nin, NeMid);

    % Inputs/stimuli -> midE

    for k = 1:Nstim

        Aconn = tempconn;
        Aconn(:, subLayerIndicesInM(1, k):subLayerIndicesInM(2, k)) =  0.74;
    
        c = length(s.connections) + 1;
        s.connections(c).direction = [IOping{k}.populations(1).name, '->', pingM.populations(1).name];
        s.connections(c).source = IOping{k}.populations(1).name;
        s.connections(c).target = pingM.populations(1).name;
        s.connections(c).mechanism_list={'input_connection'};
        s.connections(c).parameters={'input_amp',input_amp,'netcon',Aconn};
    
        c = length(s.connections) + 1;
        s.connections(c).direction = [IOping{k}.populations(2).name, '->', pingM.populations(1).name];
        s.connections(c).source = IOping{k}.populations(2).name;
        s.connections(c).target = pingM.populations(1).name;
        s.connections(c).mechanism_list={'input_connection'};
        s.connections(c).parameters={'input_amp',input_amp,'netcon',Aconn};

    end
    
    % supEmidE
    c = length(s.connections)+1;
    s.connections(c).direction = [pingS.populations(1).name, '->', pingM.populations(1).name];
    s.connections(c).source = pingS.populations(1).name;
    s.connections(c).target = pingM.populations(1).name;
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_MS,'tauAMPA',tauAMPA_SM,'netcon',KsupEmidE};

    % midEsupE
    c = length(s.connections)+1;
    s.connections(c).direction = [pingM.populations(1).name, '->', pingS.populations(1).name];
    s.connections(c).source = pingM.populations(1).name;
    s.connections(c).target = pingS.populations(1).name;
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_MS,'tauAMPA',tauAMPA_MS,'netcon',KmidEsupE};

    % midEdeepE
    c = length(s.connections)+1;
    s.connections(c).direction = [pingM.populations(1).name, '->', pingD.populations(1).name];
    s.connections(c).source = pingM.populations(1).name;
    s.connections(c).target = pingD.populations(1).name;
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_MD,'tauAMPA',tauAMPA_MD,'netcon',KmidEdeepE};
    
    % deepEmidE
    c = length(s.connections)+1;
    s.connections(c).direction = [pingD.populations(1).name, '->', pingM.populations(1).name];
    s.connections(c).source = pingD.populations(1).name;
    s.connections(c).target = pingS.populations(1).name;
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_MD,'tauAMPA',tauAMPA_DM,'netcon',KdeepEmidE};

    % Outputs: deepE is the recommended output layer.
    y = s;
    fprintf("\n->Initialization of dsModel ""%s"" is done. \n", populationName);

end