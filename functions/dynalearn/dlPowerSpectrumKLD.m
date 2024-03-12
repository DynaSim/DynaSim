function d = dlPowerSpectrumKLD(dlObj, opts)

    dlPotentialIndices = endsWith(dlObj.dlVariables, opts.id);
    dlPotentials = dlObj.dlOutputs(dlPotentialIndices);
    dlFs = floor(1000/(dlObj.dldT * dlObj.dlDownSampleFactor));  
    dlQ = opts.target;

    n = size(dlPotentials, 2);
    x = dlPotentials{1, 1};
    L = size(x, 1);

    dlLf = 1 + floor((opts.lf / dlFs)*L);
    dlHf = floor((opts.hf / dlFs)*L);
    y = zeros(n, L);
            
    for k = 1:n
        
        x = dlPotentials{1, k};
        N = size(x, 2);
        x = mean(x, 2);

        x = detrend(x);
        [Y, ~] = dlSpectrum(x, dlFs);
        y(k, :) = Y*N;
        
    end
    
    dlP = mean(y, 1);
    dlP = dlP(dlLf:dlHf);
    d = dlKullbackLeiblerDivergence(dlP, dlQ);

end

