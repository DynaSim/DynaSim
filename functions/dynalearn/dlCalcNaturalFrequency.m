function y = dlCalcNaturalFrequency(dlOutput, dldt, dlDownSampleFactor)

    L = size(dlOutput, 1);
    Lnq = floor(L/2)+1;
    N = size(dlOutput, 2);
    Fs = (1 / (dldt*dlDownSampleFactor)) * 1000;

    x = dlOutput;
    x = mean(x, 2);
    x = detrend(x);

    Y = fft(x); 
    P2 = abs(Y/L);
    P1 = P2(1:Lnq);
    f = Fs*(0:Lnq)/L;

    P1 = smooth(P1, min(N+1, ceil(L^0.25)));
    y = f(find(P1 == max(P1), 1));

    % plot(f, P1);title("FFT");xlabel("Frequency");ylabel("Power amplitude");

end
