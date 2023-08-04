function [y, l] = dlRPowerSpectrum(dlObj, opts)

    dlPotentialIndices = endsWith(dlObj.dlVariables, opts.id);
    dlPotentials = dlObj.dlOutputs(dlPotentialIndices);
    dlLabels = dlObj.dlVariables(dlPotentialIndices);

    n = size(dlPotentials, 2);
    y = zeros(n, 1);
            
    Fnq = 0.5 * (1 / (dlObj.dldT*dlObj.dlDownSampleFactor)) * 1000;
    
    for k = 1:n
        
        x = dlPotentials{1, k};
        L = size(x, 1);
        N = size(x, 2);
        x = mean(x, 2);

        x = detrend(x);
        Y = fft(x);
        Lnq = floor(L/2)+1;

        P2 = abs(Y/L);
        P1 = P2(1:Lnq);
        P1 = smooth(P1, min(N+1, ceil(L^0.25)));

        lf1 = ceil(Lnq*opts.lf1/Fnq);
        hf1 = floor(Lnq*opts.hf1/Fnq);
        lf2 = ceil(Lnq*opts.lf2/Fnq);
        hf2 = floor(Lnq*opts.hf2/Fnq);

        yf1 = mean(P1(lf1:hf1));
        yf2 = mean(P1(lf2:hf2));

        y(k) = mean(yf1)/mean(yf2);
        
    end
    
    l = dlLabels;

end

