function [y, l] = dlEPowerSpectrum(dlObj, opts)

    dlPotentialIndices = contains(dlObj.dlVariables, '_V');
    dlPotentialIndices(1) = 1;
    dlPotentials = dlObj.dlOutputs(dlPotentialIndices);
    dlLabels = dlObj.dlVariables(dlPotentialIndices);

    n = size(dlPotentials, 2);
    y = zeros(n, 1);
            
    dtf = ceil(1 / (dlObj.dldT*dlObj.dlDownSampleFactor));

    lf1 = opts.lf1*dtf;
    hf1 = opts.hf1*dtf;
    
    for k = 1:n
        
        x = dlPotentials{1, k};
        ffts = abs(fft(mean(x, 2))); 
        yf1 = smooth(ffts(lf1:hf1));

        y(k) = mean(yf1);
        
    end
    
    l = dlLabels;

end
