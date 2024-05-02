function neuralNetworkCircuit = dlCorticalColumn(specs)
    
    eqns={
      'dV/dt = (Iapp)*(cos(2*pi*omega*t))*(t > t0 & t < t1) + @current + noise*(rand(1, N_pop)); noise = 0; t0 = 1000; t1 = 1500; Iapp = 0; omega = 0; V(0) = -70*rand(1, N_pop);'
    };

    s=[];
    pcnt = 0;

    for i = 1:specs.layers

        for j = 1:specs.celltypesCount

            pcnt = pcnt + 1;
            pname = specs.celltypes(j);
            pnoise = specs.noises(j);
            psize = sepcs.counts(i, j);

            s.populations(pcnt).name = pname + num2str(i);
            s.populations(pcnt).size = psize;
            s.populations(pcnt).equations = eqns;
            s.populations(pcnt).mechanism_list = {'iNa', 'iK', 'ileak'};
            s.populations(pcnt).parameters = {'gNa', 120, 'gK', 36, 'gleak', .6, 'noise', pnoise};
        
        end

    end

    for i = 1:pcnt

        for j = 1:pcnt

            psource = s.populations(i).name;
            psync = s.populations(j).name;
            
            if contains(psource, "PV")

                s.connections(1).direction='INfast->ES';
                s.connections(1).mechanism_list={'iGABAa'};
                s.connections(1).parameters={'tauD', 11,'gGABAa', 0.2, 'netcon', 'zeros(n_pre, n_post)'};
    
            elseif contains(psource, "CB")

            elseif contains(psource, "CR")

            else

            end

        end

    end
    
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