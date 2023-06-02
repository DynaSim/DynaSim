function y = dlCalcNaturalFrequency(dlOutput, dldt, dlDownSampleFactor)

    L = size(dlOutput, 1);
    Fs = (1 / (dldt*dlDownSampleFactor)) * 1000;

    x = dlOutput;
    x = mean(x, 2);
    x = detrend(x);

    Y = fft(x); 
    P2 = abs(Y/L);
    P1 = P2(1:floor(L/2)+1);
    f = Fs*(0:floor(L/2))/L;

    y = f(find(P1 == max(P1), 1));

end
