function dlSpectrogramPlot(X, timeW, freqW, overlap, fmax, fs)

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

            m = size(X, 1);
            tmax = size(X, 2);

            tW = (timeW - overlap)*(fs/1000);
            tWx = (timeW)*(fs/1000);
            tB = floor(tmax / tW)-1;
            fB = floor(fmax / freqW);

            y = zeros(m, tB, fB);

            for i = 1:m

                for j = 1:tB

                    tK = max(j*tW - tWx, 1):min(j*tW, tmax);
                    tempX = conv(X(i, tK), exp(linspace(1, 0, 7)), "same");
                    tempF = dlSpectrum(tempX, fs, fmax, fB);
                    y(i, j, :) = tempF;

                end

            end

            t = linspace(0, tmax, tB);
            f = linspace(0, fmax, fB);
            figure('Position', [0, 0, 1700, 1400]);

            sG = squeeze(mean(y, 1));
            subplot(1, 1, 1);
            imagesc(sG', "XData", t, "YData", f);
            xlabel("Time (ms)");ylabel("Freq (Hz)");
            colormap("jet");


            sgtitle("Spectrogram");

        end