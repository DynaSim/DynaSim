function y = dlSpectrogramPlot(X, timeW, freqW, overlap, fmax, fs, Skernel)

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

            tK = max(j*tW - tWx, 1):min(j*tW, tmax);
            tempT = exp(-(linspace(-.5, 4.5, kernelSize).^2));
            tempX = conv(X(i, tK), tempT/sum(tempT), "same");
            tempF = dlSpectrum(tempX, fs, fmax, fB);
            y(i, j, :) = tempF;

        end

    end

    t = linspace(0, tmax, tB);
    f = linspace(0, fmax, fB);
    figure('Position', [0, 0, 1700, 1400]);

    sG = squeeze(mean(y, 1));
    y = sG';
    subplot(1, 1, 1);
    imagesc(y, "XData", t, "YData", f);
    
    xlabel("Time (ms)");
    ylabel("Freq (Hz)");
    colormap("jet");
    sgtitle("Spectrogram");

end