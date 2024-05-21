function d = dlPowerSpectrumLogCorrelation2(dlObj, opts)

    dlFs = floor(1000/(dlObj.dldT * dlObj.dlDownSampleFactor));  
    dlQ = opts.target;

    x = dlObj.dlSignals;
    dlP = dlSpectrogramPlot(x, 200, 1, 180, opts.hf, dlFs, 5, 0);

    d = dlLogCorrelationDivergence(dlP, dlQ);

end

