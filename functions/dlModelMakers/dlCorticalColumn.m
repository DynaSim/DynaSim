function neuralNetworkCircuit = dlCorticalColumn(name, Specs)

    eqns={
      'dV/dt = (Iapp)*(cos(2*pi*omega*t))*(t > t0 & t < t1) + @current + noise*(rand(1, N_pop)); noise = 0; t0 = 1000; t1 = 1500; Iapp = 0; omega = 0; V(0) = -70*rand(1, N_pop);'
    };
    
    nES = 70;
    nIF = 10;
    nIS = 10;
    nIV = 10;
    
    gGABAslow = 0.25; % Receptor conductances, will be modulated by synaptic conn.
    gGABAfast = 0.2;
    gAMPA = 0.1;
    
    tauGABAslow = 100*(1 + .1*randn(1)); % Receptor synaptic time constants; SST->~10Hz ~(1/2t = 1/100)
    tauGABAfast = 10*(1 + .1*randn(1)); % PV->~100Hz (1/10)
    tauAMPAnear = 20*(1 + .1*randn(1)); % PY->~50.0Hz (1/20)
    tauAMPAfar = 50*(1 + .1*randn(1)); % PY->~25.0Hz (1/40)
    
    eNoise = 25;
    fNoise = 25;
    sNoise = 25;
    
    cCoeff = .01;
    
    kEE = ones(nES, nES) * cCoeff;
    kEIf = ones(nES, nIF) * cCoeff;
    kEIs = ones(nES, nIS) * cCoeff;
    
    kIfE = ones(nIF, nES) * cCoeff;
    kIfIf = ones(nIF, nIF) * cCoeff;
    kIfIs = ones(nIF, nIS) * cCoeff;
    
    kIsE = ones(nIS, nES) * cCoeff;
    kIsIf = ones(nIS, nIF) * cCoeff;
    kIsIs = ones(nIS, nIS) * cCoeff;
    
    s=[];
    s.populations(1).name='ES';
    s.populations(1).size=nES;
    s.populations(1).equations=eqns;
    s.populations(1).mechanism_list={'iNa','iK', 'ileak'};
    s.populations(1).parameters={'gNa',120,'gK',36, 'gleak', .6,'noise',eNoise};
    
    s.populations(2).name='INfast';
    s.populations(2).size=nIF;
    s.populations(2).equations=eqns;
    s.populations(2).mechanism_list={'iNa','iK', 'ileak'};
    s.populations(2).parameters={'gNa',120,'gK',36, 'gleak', .5,'noise',fNoise};
    
    s.populations(3).name='INslow';
    s.populations(3).size=nIS;
    s.populations(3).equations=eqns;
    s.populations(3).mechanism_list={'iNa','iK', 'ileak'};
    s.populations(3).parameters={'gNa',120,'gK',36, 'gleak', .5,'noise',sNoise};
    
    s.connections(1).direction='INfast->ES';
    s.connections(1).mechanism_list={'iGABAa'};
    s.connections(1).parameters={'tauD',tauGABAfast,'gGABAa',gGABAfast, 'netcon', kIfE};
    
    s.connections(2).direction='ES->INfast';
    s.connections(2).mechanism_list={'iAMPA'};
    s.connections(2).parameters={'tauD',tauAMPAnear,'gAMPA',gAMPA, 'netcon', kEIf};
    
    s.connections(3).direction='INslow->ES';
    s.connections(3).mechanism_list={'iGABAa'};
    s.connections(3).parameters={'tauD',tauGABAslow,'gGABAa',gGABAslow, 'netcon', kIsE};
    
    s.connections(4).direction='ES->INslow';
    s.connections(4).mechanism_list={'iAMPA'};
    s.connections(4).parameters={'tauD',tauAMPAfar,'gAMPA',gAMPA, 'netcon', kEIs};
    
    s.connections(5).direction='INslow->INfast';
    s.connections(5).mechanism_list={'iGABAa'};
    s.connections(5).parameters={'tauD',tauGABAslow,'gGABAa',gGABAslow, 'netcon', kIsIf};
    
    s.connections(6).direction='INslow->INslow';
    s.connections(6).mechanism_list={'iGABAa'};
    s.connections(6).parameters={'tauD',tauGABAslow,'gGABAa',gGABAslow, 'netcon', kIsIs};
    
    s.connections(7).direction='INfast->INslow';
    s.connections(7).mechanism_list={'iGABAa'};
    s.connections(7).parameters={'tauD',tauGABAfast,'gGABAa',gGABAfast, 'netcon', kIfIs};
    
    s.connections(8).direction='INfast->INfast';
    s.connections(8).mechanism_list={'iGABAa'};
    s.connections(8).parameters={'tauD',tauGABAfast,'gGABAa',gGABAfast, 'netcon', kIfIf};
    
    s.connections(9).direction='ES->ES';
    s.connections(9).mechanism_list={'iAMPA'};
    s.connections(9).parameters={'tauD',tauAMPAnear,'gAMPA',gAMPA, 'netcon', kEE};

    neuralNetworkCircuit = s;

end