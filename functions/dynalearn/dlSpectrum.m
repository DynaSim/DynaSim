function [y, f] = dlSpectrum(x, fs, fmax, fcnt, smoothing)

    N = length(x);

    if ~exist('fs', 'var')

        fs = 1000;
        disp("No sampling rate was specificed. Default Fs = 1000Hz ");

    end

    if ~exist('fmax', 'var')

        fmax = fs/2;

    end

    if ~exist('fcnt', 'var')

        fcnt = N;

    end

    if ~exist('smoothing', 'var')

        smoothing = 1;

    end

    t1 = linspace(0, fs/2, ceil(N/2)-1);
    t2 = linspace(0, fmax, fcnt);

    if smoothing

        x = smooth(x, ceil(fs / (fmax)));

    end

    p1 = fft(x);
    p2 = abs(p1(2:ceil(N/2)));
    y = interp1(t1, p2, t2);
    f = t2;

end