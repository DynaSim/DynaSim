function y = dlSinglePopulationNS(ModelParameters, populationName)

    try

        Ne = ModelParameters.Size;
        NoiseRate = ModelParameters.Noise;

    catch

        fprintf("->Failed to make model. \n->Inputs must be a ModelParameters struct with Size and Noise specified, plus populationName as a string");
        error("->Invalid use of the function.");

    end

    totalNeuronCount = Ne;

    fprintf("\n->Initialization of dlSinglePopulationNS Model:\n-->Total neuron count = %d ", totalNeuronCount);
    fprintf("\n-->Overall noise rate = %.4f", NoiseRate); % Randomness / Stochasticity
    fprintf("\n->Population name is %s", populationName); % Name tag or suffix for all names of this dsModel

    k1 = 0.2; % Diff. for normal weights (uniform random)

    ModelName = populationName;
    populationName = char("_" + populationName);

    % Connectivity matrices
    gBase = 7/(Ne^0.5);

    % neuronal dynamics
    eqns = 'dV/dt = (Iapp + @current + noise*randn(1, Npop))/C; Iapp=0; noise=0; C=1; V(0) = -65 + rand(1, Npop)*0;';
    % eqns_apical = 'dV/dt = (Iapp + Icable + @current + noise*randn(1, Npop))/C; Iapp=0; noise=0; C=1; V(0) = -65 + rand(1, Npop)*0;';
    %eqns_stimuli = 'dinp/dt = (rand(1) + 4.5)*(20*(exp(- (t - t1).^2) - exp(- (t - t2).^2)) + noise*randn(1, Npop))/C; f1=4; t1=10; t2=100; noise=0; C=1; inp(0) = -60 - rand(1, Npop)*0;';
    % eqns_stimuli = 'dinp/dt = (noise*rand(1) + 1)*(1*(exp(- (t - t1).^2) - exp(- (t - t2).^2)) + noise*randn(1, Npop)); t1=100; t2=200; noise=0; inp(0) = rand(1, Npop)*0;';

    % cell type
    ctx_cells = {'iNa','iK', 'iCa', 'iL'};
    cell_type = ctx_cells;

    % BUILDING
    % Structures: SupERFICIAL LAYER (1-2-3)
    dModel=[];

    % E-cells-PYR
    dModel.populations(1).name = ['ESOMA', populationName];
    dModel.populations(1).size = Ne;
    dModel.populations(1).equations = eqns;
    dModel.populations(1).mechanism_list = cell_type;
    dModel.populations(1).parameters = {'Iapp', 0,'noise', NoiseRate*2};

    % I-cells-PV
    dModel.populations(2).name = ['EDEND', populationName];
    dModel.populations(2).size = Ne;
    dModel.populations(2).equations = eqns;
    dModel.populations(2).mechanism_list = cell_type;
    dModel.populations(2).parameters = {'Iapp', 0,'noise', NoiseRate};

    dModel.connections(1).direction = [dModel.populations(1).name, '->', dModel.populations(2).name];
    dModel.connections(1).source = dModel.populations(1).name;
    dModel.connections(1).target = dModel.populations(2).name;
    dModel.connections(1).mechanism_list = {'iCOM'};
    dModel.connections(1).parameters = {'gCOM',gBase};

    dModel.connections(2).direction = [dModel.populations(2).name, '->', dModel.populations(1).name];
    dModel.connections(2).source = dModel.populations(2).name;
    dModel.connections(2).target = dModel.populations(1).name;
    dModel.connections(2).mechanism_list = {'iCOM'};
    dModel.connections(2).parameters = {'gCOM',gBase};

    y = dModel;
    fprintf("\n->Initialization of dsModel ""%s"" is done. \n", ModelName);

end