%% Load model logs

clear;clc;

ratioLogs = cell(50, 1);
lfps = zeros(50, 750, 24290);
responses = cell(50, 1);

for i = 1:1

    m = DynaLearn(); % ~ 1sec
    m = m.dlLoad(char("dlModels/dlTCTModels" + string(i))); % ~ 10sec, New larger model; keeping track of its activity in Gamma/Beta **
   
    ratioLog1 = {m.dlCustomLog(1, :)};
    ratioLog1 = cell2mat(ratioLog1{1});
    ratioLogs_1{i} = ratioLog1';

    ratioLog2 = {m.dlCustomLog(2, :)};
    ratioLog2 = cell2mat(ratioLog2{1});
    ratioLogs_2{i} = ratioLog2';

    ratioLog3 = {m.dlCustomLog(3, :)};
    ratioLog3 = cell2mat(ratioLog3{1});
    ratioLogs_3{i} = ratioLog3';

    ratioLog4 = {m.dlCustomLog(4, :)};
    ratioLog4 = cell2mat(ratioLog4{1});
    ratioLogs_4{i} = ratioLog4';

    lfpcl = {m.dlCustomLog(5, :)};
    lfps(i, :, :) = cell2mat(lfpcl{1});

    rsl = {m.dlCustomLog(6, :)};
    responses{i} = cell2mat(rsl{1});

end

%% To matrix

clc;

ratioLogs1 = cell2mat(ratioLogs_1);
ratioLogs2 = cell2mat(ratioLogs_2);
ratioLogs3 = cell2mat(ratioLogs_3);
ratioLogs4 = cell2mat(ratioLogs_4);

ratioLogs1 = ratioLogs1'; % Theta
ratioLogs2 = ratioLogs2'; % Alpha
ratioLogs3 = ratioLogs3'; % Beta
ratioLogs4 = ratioLogs4'; % Gamma

save("aa.mat", 'm', 'ratioLogs1', 'ratioLogs2', 'ratioLogs3', 'ratioLogs4', 'responses');

%% Averages

clc;clear;
load("aa.mat");

modelType = "Inhibitory Excluded";
bandLabel = "Gamma";
qq = ratioLogs4;

ratioLogsPooled = 100*(qq - mean(qq, 'all')) / mean(qq, 'all');
ratioLogsAvg = zeros(size(ratioLogs1));

for i = 1:450

    ratioLogsPooled(i, :) = conv(ratioLogsPooled(i, :), [1 1 1]/3, "same");
%     ratioLogsPooled(i, :) = ratioLogsPooled(i, :) - mean(ratioLogsPooled(i, :));

end

for i = 1:45

    q = ratioLogsPooled(i:10:90, :);
    ratioLogsAvg(i, :) = mean(q, 1);
%     ratioLogsAvg(i, :) = ratioLogsAvg(i, :) - ratioLogsAvg(i, 5);

end

for i = 1:50

    q = responses{i};
    responses{i} = q(1:694);
    
end

%% Ratio change (relative) in trials during predictive task training

clc;

w = 100;
cnt = 1;
tlab = ["U0", "A", "U1", "B", "U2", "C", "U3", "A", "U4", "B", "U5", "C", "U6"];
figure("Position", [0 0 1500 1000]);

set1 = 2:9;
% ratioLogPlot = ratioLogsAvg(:, 5:590);

for i = 1:8

    for l = 0:6

        subplot(2, 4, i);
        fill([l*w, l*w, l*w+w, l*w+w], [-1e+5, 1e+5, 1e+5, -1e+5], [sin(l*0.2), 1, cos(l*0.2)], 'HandleVisibility','off');
        hold('on');

    end

end

dispLabels = m.dlCustomLogLabel{1};

for i = set1

    subplot(2, 4, mod(cnt-1, 8)+1);
    legend("Location", "southwest");

    x = 5:694;
    q = ratioLogsPooled(i:10:90, x);
                                          % Create Independent Variable
    y = q;                                  % Create Dependent Variable 6Experimentsz Data
    N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
    yMean = mean(y, 1);                                    % Mean Of All Experiments At Each Value Of ,x 
    ySEM = std(y)/sqrt(N);                              % Compute 2Standard Error Of The Mean* Of All Experiments At Each Value Of .x1
    y1 = yMean + ySEM;
    y2 = yMean - ySEM;
%     CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
%     yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

    ql = min([y1, y2]);
    qu = max([y1, y2]);

    for l = 0:6
        text((l+0.5)*w, qu+.1, tlab(l+1));
    end

    ylim([ql-.2 qu+.2]);
    title(dispLabels{i});
    
    plot(x, yMean, '-', 'DisplayName', dispLabels{i});                                 % Plot Mean Of All Experiments
    hold("on");
%     plot(x, yCI95+yMean,'-r')                           % Plot 95% Confidence Intervals Of All Experiments
    patch([x, fliplr(x)], [y1, fliplr(y2)], 'b', 'EdgeColor','g', 'FaceAlpha',0.25, 'DisplayName', 'CI-95%');

    xlim([min(x), max(x)]);
    hold("on");
    grid("on");
    cnt = cnt + 1;

end

sgtitle(bandLabel + " power change (Percentage) across trials (Average across 10 different " + modelType + " models)");

% Ratio change (relative) in trials during predictive task training

clc;

w = 100;
cnt = 1;
tlab = ["P", "U", "P", "U", "P", "U"];
f = figure("Position", [0 0 1500 1000]);

set1 = 2:9;

ratioLogPlot = ratioLogsAvg(:, 5:590);

for i = 1:8

    for l = 1:3

        subplot(2, 4, i);
        fill([l*w, l*w, l*w+w, l*w+w], [-1e3, 1e3, 1e3, -1e3], [sin(l*0.45), 1, cos(l*0.45)], 'HandleVisibility','off');
        hold('on');

    end

end

dispLabels = m.dlCustomLogLabel{1};

for i = set1

    subplot(2, 4, mod(cnt-1, 8)+1);
%     plot(ratioLogPlot(i, :), 'DisplayName', dispLabels{i});
    legend("Location", "southwest");

%     for l = 0:5
%         text((l+0.5)*w, qu+0.5, tlab(l+1));
%     end

    x = 101:300;
    q = ratioLogsPooled(i:10:90, x);
                                          % Create Independent Variable
    y = q;                                  % Create Dependent Variable 6Experimentsz Data
    N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
    yMean = mean(y, 1);                                    % Mean Of All Experiments At Each Value Of ,x 
    ySEM = std(y)/sqrt(N);                              % Compute 2Standard Error Of The Mean* Of All Experiments At Each Value Of sxq
    y1 = yMean + ySEM;
    y2 = yMean - ySEM;
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

    ql = min([y1, y2]);
    qu = max([y1, y2]);

    ylim([ql-.2 qu+.2]);
    xlim([min(x), max(x)])
    title(dispLabels{i});
    
    plot(x, yMean, '-', 'DisplayName', dispLabels{i});                                  % Plot Mean Of All Experiments
    hold("on");
%     plot(x, yCI95+yMean,'-r')                           % Plot 95% Confidence Intervals Of All Experiments
    patch([x, fliplr(x)], [y1, fliplr(y2)], 'b', 'EdgeColor','g', 'FaceAlpha',0.25, 'DisplayName', 'CI-95%');

    hold("on");
    grid("on");
    cnt = cnt + 1;

end

sgtitle(bandLabel + " power band change % (Average across 10 different " + modelType + " models/ Predictable block (green) to unpredictable (yellow) transitions)");

%%

clc;
Currentsize = 30;noise_rate=5.0;

for i = 1

    %%% Create model parameters struct
    ModelParametersPFC = struct();
    
    %%% Area PFC layer sizes (relative)
    ModelParametersPFC.NeSuperficial = ceil(0.3*Currentsize);
    ModelParametersPFC.NSomSuperficial = ceil(0.03*Currentsize);
    ModelParametersPFC.NPvSuperficial = ceil(0.03*Currentsize);
    ModelParametersPFC.NeMid = ceil(0.14*Currentsize);
    
    ModelParametersPFC.NSomMid = 0;
    
    ModelParametersPFC.NPvMid = ceil(0.03*Currentsize);
    ModelParametersPFC.NeDeep = ceil(0.45*Currentsize);
    ModelParametersPFC.NSomDeep = ceil(0.04*Currentsize);
    ModelParametersPFC.NPvDeep = ceil(0.01*Currentsize);
    
    ModelParametersPFC.Nin = 6;
    ModelParametersPFC.Nout = 6;
    ModelParametersPFC.NoiseRate = noise_rate;
    ModelParametersPFC.Nstim = 3;
    ychlabels = ["Time"];
    
    for i = 1:ModelParametersPFC.NeSuperficial
        ychlabels = [ychlabels, "supE"];
    end
    for i = 1:ModelParametersPFC.NSomSuperficial
        ychlabels = [ychlabels, "supISOM"];
    end
    for i = 1:ModelParametersPFC.NPvSuperficial
        ychlabels = [ychlabels, "supIPV"];
    end
    
    for i = 1:ModelParametersPFC.NeMid
        ychlabels = [ychlabels, "midE"];
    end
    for i = 1:ModelParametersPFC.NPvMid
        ychlabels = [ychlabels, "midIPV"];
    end
    
    for i = 1:ModelParametersPFC.NeDeep
        ychlabels = [ychlabels, "deepE"];
    end
    for i = 1:ModelParametersPFC.NSomDeep
        ychlabels = [ychlabels, "deepISOM"];
    end
    for i = 1:ModelParametersPFC.NPvDeep
        ychlabels = [ychlabels, "deepIPV"];
    end
    
    ychlabels = flip(ychlabels);

end

popsize = size(ychlabels, 2);
trial = 202;
k = 9;

lfptemp = {m.dlCustomLog{5, :}};
lfptemp = cell2mat(lfptemp);

%

figure("Position", [0 0 1500 1000]);
subplot(1, 1, 1);

for i = 1:popsize

%     y = lfps(k, :, i + popsize*trial);
    y = lfptemp(:, i + popsize*trial);
    y = y - min(y);
    y = y / max(y);
    offsetLFPplot = (popsize+1-i) * 1;

    plot(y + offsetLFPplot);hold("on");
    
end

% set(gca, 'YDir','reverse');
yticklabels(ychlabels);
yticks(1:popsize);
ylim([-1 popsize+1])

%% Accuracy

clc;
a = cell2mat(responses);
b = a;

b(:, 1:100) = (b(:, 1:100) == 1);
b(:, 101:200) = (b(:, 101:200) ~= 1);

b(:, 201:300) = (b(:, 201:300) == 2);
b(:, 301:400) = (b(:, 301:400) ~= 2);

b(:, 401:500) = (b(:, 401:500) == 3);
b(:, 501:end) = (b(:, 501:end) ~= 3);

for i = 1:size(b, 1)

    b(i, :) = dlMpool1d(b(i, :), 10);

end

b1 = b(11:20, :);
b2 = b(41:50, :);
c1 = mean(b1, 1);
c2 = mean(b2, 1);

figure("Position", [0 0 1500 1000]);
subplot(2, 1, 1);

for l = 0:6

    subplot(2, 1, 1);
    fill([l*w, l*w, l*w+w, l*w+w], [-1e3, 1e3, 1e3, -1e3], [sin(l*0.25), 1, cos(l*0.25)], 'HandleVisibility','off');
    hold('on');

end

plot(c1, 'DisplayName', 'Normal models');hold('on');
plot(c2, 'DisplayName', 'Excluded inhibitory layers from learning');
title("Avg. trial accuracy across 600 trials for 10 " + modelType + " models.");
grid("on");ylim([0 1]);ylabel("Accuracy");
legend('Location', 'southeast');
subplot(2, 1, 2);imagesc(b1);title("Response of model (yellow=Correct, blue=Incorrect)");
xlabel("Trial");ylabel("Model no.");
yticklabels(1:20);
% yticks(1:10);


%% End
