function y = dlSingleLayer(ModelParameters, populationName)

    NE = ModelParameters.NE;
    NSom = ModelParameters.NSom;
    NPv = ModelParameters.NPv;

    Nstim = ModelParameters.Nstim;
    NoiseRate = ModelParameters.NoiseRate;

    fprintf("\n>Initialization of dlSingleLayer Model: ");
    fprintf("\n--> Excitatory neurons count = %d , SOM inhibitory = %d , PV inhibitory = %d ", Ne, NSom, NPv);

    fprintf("\n-->Input connections count (terminal) size = %d ", NStim); % Inputs / Stimuli
    fprintf("\n-->Overall noise rate = %.4f", NoiseRate); % Randomness / Stochasticity
    fprintf("\n--->Population name is %s", populationName); % Name tag or suffix for all names of this dsModel

    k1 = 0.15; % Diff. for normal weights (uniform random)
    k2 = 0.45; % Min connectivity weight

    NeAvg = Ne + NSom + NPv;

    populationName = ['x', populationName];
    % Connectivity matrices

    % sE->sIsom
    KsupEsupSom = k1*rand(NE, NSom) + 0.7;
    % sE->sIpv
    KsupEsupPv = k1*rand(NE, NPv) + 0.7;

    % sIsom->sE
    KsupSomsupE = k1*rand(NSom, NE) + k2;
    % sIsom->sIpv
    KsupSomsupPv = k1*rand(NSom, NPv) + k2;

    % sIpv->sE
    KsupPvsupE = k1*rand(NPv, Ne) + 0.7;
    % sIpv->sIsom
    KsupPvsupSom = k1*rand(NPv, NSom) + k2;
    
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
    tauGABA_Som = 151; % ms, decay time constant of inhibition for beta (around 25Hz)
    tauGABA_PvSup = 57;

    tauAMPA_E = 171;
    
    gBase = 7/(NeAvg^0.5);

    % Maximal synaptic strengths
    gGABA_SE = .017*gBase; % (SOM)->E
    gGABA_PE = .011*gBase; % (PV)->E
    gGABA_II = .011*gBase; % I->I within layer


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

    s = dsCombineSpecifications(pingS, IOping{1}, IOping{2}, IOping{3});

    % connect the layers and inputs
    fprintf("\n--->Connecting separate layers and inputs:");
    tempconn = zeros(Nin, NeMid);

    % Inputs/stimuli -> E

    for k = 1:Nstim

        Aconn = tempconn;
        Aconn(:, subLayerIndicesInM(1, k):subLayerIndicesInM(2, k)) =  0.74;
    
        c = length(s.connections) + 1;
        s.connections(c).direction = [IOping{k}.populations(1).name, '->', pingS.populations(1).name];
        s.connections(c).source = IOping{k}.populations(1).name;
        s.connections(c).target = pingS.populations(1).name;
        s.connections(c).mechanism_list={'input_connection'};
        s.connections(c).parameters={'input_amp',input_amp,'netcon',Aconn};
    
        c = length(s.connections) + 1;
        s.connections(c).direction = [IOping{k}.populations(2).name, '->', pingS.populations(1).name];
        s.connections(c).source = IOping{k}.populations(2).name;
        s.connections(c).target = pingS.populations(1).name;
        s.connections(c).mechanism_list={'input_connection'};
        s.connections(c).parameters={'input_amp',input_amp,'netcon',Aconn};

    end
    
    % Outputs: deepE is the recommended output layer.
    y = s;
    fprintf("\n->Initialization of dsModel ""%s"" is done. \n", populationName);

end