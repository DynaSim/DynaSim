function [y, f] = dlSpectrum(x, fs)

    if ~exist('fs', 'var')

        fs = 1000;
        disp("No sampling rate was specificed. Default Fs = 1000Hz ");

    end

    N = length(x);
    t1 = linspace(0, 1, ceil(N/2));
    t2 = linspace(0, 1, N);

    p1 = fft(x);
    p2 = abs(p1(1:ceil(N/2)));
    p3 = interp1(t1, p2, t2);

    y = smooth(p3, ceil(N^0.5));
    y = y / sum(y);
    f = linspace(0, fs, N);

end