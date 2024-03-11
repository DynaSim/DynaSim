function [y, f] = dlSpectrum(x, fs)

    if ~exist('fs', 'var')

        fs = 1000;
        disp("No sampling rate was specificed. Default Fs = 1000Hz ");

    end

    N = length(x);
    x = detrend(x);
    t1 = linspace(0, fs/2, ceil(N/2));
    t2 = linspace(0, fs/2, N);

    p1 = fft(x);
    p2 = abs(p1(1:ceil(N/2)));
    p2 = smooth(p2, min(10, ceil(N^0.5)));
    y = interp1(t1, p2, t2);

    f = linspace(0, fs, N);

end