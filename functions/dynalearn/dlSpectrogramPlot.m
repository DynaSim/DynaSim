function [y, t, f] = dlSpectrogramPlot(X, timeW, freqW, overlap, fmax, fs, Skernel, plotFlag)

    if ~exist('timeW', 'var')

        timeW = 100;

    end

    if ~exist('freqW', 'var')

        freqW = 1;

    end

    if ~exist('overlap', 'var')

        overlap = 90;

    end

    if ~exist('fmax', 'var')

        fmax = 100;

    end

    if ~exist('fs', 'var')

        fs = 1000;

    end

    if ~exist('Skernel', 'var')

        Skernel = 10;

    end

    if ~exist('plotFlag', 'var')

        plotFlag = 1;

    end

    kernelSize = ceil(fs / 1000)*Skernel;
    tempT = exp(-abs(linspace(-.5, 2.5, kernelSize).^2));
    X = conv2(X, tempT, "same");
    X = mean(X, 1);

    [sG, ~, ~] = pspectrum(X, fs, "spectrogram", "FrequencyLimits", [0 fmax], "TimeResolution", 0.4, "OverlapPercent", 95);

    sG = sG / max(max(sG));
    y = sG';

    if plotFlag

        tmax = size(X, 2);
        tW = (timeW - overlap)*(fs/1000);
        tB = floor(tmax / tW);
        fB = floor(fmax / freqW);

        t = linspace(0, tmax, tB);
        f = linspace(0, fmax, fB);

        figure('Position', [0, 0, 1700, 1400]);
        subplot(1, 1, 1);
        imagesc(y, "XData", t, "YData", f);
        
        xlabel("Time (ms)");
        ylabel("Freq (Hz)");
        colormap("jet");
        sgtitle("Spectrogram");

    end

end