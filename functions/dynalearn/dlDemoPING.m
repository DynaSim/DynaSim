function y = dlDemoPING(Ne, Ni, Nio, noise_rate)

    fprintf("Initialization...\n");

    % Population sizes

    k1 = 0.1; % Diff. for normal weights (uniform random)
    k2 = 0.2; % Min connectivity weight
    
    % Connectivity matrices

    % E->I
    Kei = k1*rand(Ne, Ni) + 0.3; % all-to-all, connectivity from E cells to I cells; mid, sup, deep
    % I->E
    Kie = k1*rand(Ni, Ne) + 0.3; % all-to-all

    % E->E
    Kee = k1*rand(Ne, Ne) + k2; % recurrent E-to-E: mid, sup, deep
    Kii = k1*rand(Ni, Ni) + k2; % recurrent I-to-I: mid, sup, deep

    kzio = zeros(Nio, Nio);

    a1 = 1;a2 = ceil(1*Ne/2);
    b1 = ceil(1 + 1*Ne/2);b2 = ceil(2*Ne/2);

    % Time constants
    tauGABA_gamma = 4.8; % ms, decay time constant of inhibition for gamma (50Hz)
    tauAMPA = 4.8; % ms, decay time constant of fast excitation (AMPA)

    % Maximal synaptic strengths
    gAMPA_ei = .2*(20/Ne); % E->I within layer
    gAMPA_in = .2*(20/Ne);

    gAMPA_ee = 0.11*(20/Ne); % E->E within layer
    gGABAa_ie = 4*(20/Ne); % I->E within layer
    gGABAa_ii = 0.11*(20/Ne); % I->I within layer
%     noise_rate = 12;

    % neuronal dynamics
    eqns = 'dV/dt = (Iapp + @current + noise*randn(1, Npop))/C; Iapp=0; noise=0; C=1; V(0) = -rand(1, Npop)*20;';
    g_poisson = 6.4e-4;
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
    IOping.connections(2).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',kzio};

    IOping.connections(3).direction = 'I->E';
    IOping.connections(3).mechanism_list = {'iAMPActx'};
    IOping.connections(3).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',kzio};

    IOping.connections(4).direction = 'I->I';
    IOping.connections(4).mechanism_list = {'iAMPActx'};
    IOping.connections(4).parameters = {'gAMPA',gAMPA_ei,'tauAMPA',tauAMPA,'netcon',kzio};

    % create independent layers
    p = dsApplyModifications(ping,{'E','name','E1'; 'I','name','I1'}); % PING
    stim = dsApplyModifications(IOping,{'E','name','S1'; 'I','name','S2'}); % I/O layer (stimuli)

    % create full cortical specification
    s = dsCombineSpecifications(p, stim);

    % connect the layers and inputs
    fprintf("Connecting separate layers and inputs...\n");
    tempconn = zeros(Nio, Ne);
    
    % Input S1 -> E1 [1]
    Aconn = tempconn;
    Aconn(:, a1:a2) =  1;

    c = length(s.connections) + 1;
    s.connections(c).direction = 'S1->E1';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Aconn};
    
    % Input S2 -> E2 [2]
    Bconn = tempconn;
    Bconn(:, b1:b2) =  1;
    
    c = length(s.connections)+1;
    s.connections(c).direction = 'S2->E1';
    s.connections(c).mechanism_list={'iAMPActx'};
    s.connections(c).parameters={'gAMPA',gAMPA_in,'tauAMPA',tauAMPA,'netcon',Bconn};
   
    % Outputs: I [1-Ne/2] as O1
    % I [Ne/2+1-Ne] as O2
    y = s;
    
    fprintf("Initialization done.\n");

end
