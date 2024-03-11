function d = dlPowerSpectrumKLD(dlObj, opts)

    dlPotentialIndices = endsWith(dlObj.dlVariables, opts.id);
    dlPotentials = dlObj.dlOutputs(dlPotentialIndices);
    dlQ = opts.target;
    dlFs = floor(1000/(dlObj.dldT * dlObj.dlDownSampleFactor));    

    n = size(dlPotentials, 2);
    x = dlPotentials{1, 1};
    L = size(x, 1);
    Lnq = floor(L/2);

    dlLf = 1 + floor((opts.lf / dlFs)*L);
    dlHf = floor((opts.hf / dlFs)*L);
    y = zeros(n, Lnq);
            
    for k = 1:n
        
        x = dlPotentials{1, k};
        N = size(x, 2);
        x = mean(x, 2);

        x = detrend(x);
        Y = fft(x);

        P2 = abs(Y/L);
        P1 = P2(1:Lnq);
        P1 = smooth(P1, min(N, ceil(L^0.5)));

        y(k, :) = P1*N;
        
    end
    
    dlP = mean(y, 1);
    dlP = dlP(dlLf:dlHf);
    d = dlKullbackLeiblerDivergence(dlP, dlQ);

end

