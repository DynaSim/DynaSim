function d = dlPowerSpectrumLogCorrelation(dlObj, opts)

    dlFs = floor(1000/(dlObj.dldT * dlObj.dlDownSampleFactor));  
    dlQ = opts.target;

    x = dlObj.dlSignals;
    n = size(x, 1);
    L = length(dlQ);

    % dlLf = 1 + floor((opts.lf / dlFs)*L);
    y = zeros(n, L);
            
    for k = 1:n
        
        Y = dlSpectrum(x(k, :), dlFs, opts.hf, L);
        y(k, :) = Y;
        
    end
    
    dlP = mean(y, 1);
    d = dlLogCorrelationDivergence(dlP, dlQ);

end

