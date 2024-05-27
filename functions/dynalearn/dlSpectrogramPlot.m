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

        Skernel = 5;

    end

    if ~exist('plotFlag', 'var')

        plotFlag = 1;

    end

    m = size(X, 1);
    tmax = size(X, 2);

    tW = (timeW - overlap)*(fs/1000);
    tWx = (timeW)*(fs/1000);
    tB = floor(tmax / tW);
    fB = floor(fmax / freqW);

    kernelSize = ceil(fs / 1000)*Skernel;
    y = zeros(m, tB, fB);

    for i = 1:m

        for j = 1:tB

            lt = max(j*tW - tWx, 1);
            rt = max(j*tW, 1);
            tK = lt:rt;
            tempX = X(i, tK);
            tempT = exp(-abs(linspace(-.5, 2.5, kernelSize).^2));
            tempX = conv(tempX, tempT/sum(tempT), "same");
            tempF = dlSpectrum(tempX, fs, fmax, fB, 0);
            y(i, j, :) = tempF;

        end

    end

    t = linspace(0, tmax, tB);
    f = linspace(0, fmax, fB);
    sG = squeeze(mean(y, 1));
    y = sG';

    if plotFlag

        figure('Position', [0, 0, 1700, 1400]);
        subplot(1, 1, 1);
        imagesc(y, "XData", t, "YData", f);
        
        xlabel("Time (ms)");
        ylabel("Freq (Hz)");
        colormap("jet");
        sgtitle("Spectrogram");

    end

end