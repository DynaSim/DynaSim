function d = dlPowerSpectrumLogCorrelation2(dlObj, opts)

    dlFs = floor(1000/(dlObj.dldT * dlObj.dlDownSampleFactor));  
    dlQ = opts.target;
    x = (dlObj.dlSignals > -10);
    dlP = dlSpectrogramPlot(x, 400, 1, 380, opts.hf, dlFs, 10, 0);

    d = dlLogCorrelationDivergence2(dlP, dlQ);

end

