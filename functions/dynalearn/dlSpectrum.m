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

    if fmax*2 > fs

        fmax = fs/2;

    end

    t2 = linspace(0, fmax, fcnt);

    if smoothing

        x = smooth(x, 2*ceil((N / (fmax)).^0.5));

    end


    [p1, t1, ~] = pspectrum(x, fs, "spectrogram", "FrequencyLimits", [0 fmax], "TimeResolution", 0.4, "OverlapPercent", 95);
    p2 = mean(p1, 2);
    y = interp1(t1, p2, t2);
    f = t2;

end