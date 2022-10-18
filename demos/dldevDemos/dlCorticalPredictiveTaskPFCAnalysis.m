%% Load model logs

% clear;clc;

ratioLogs = cell(10, 1);
lfps = zeros(10, 900, 48790);

for i = 1:7

    m = DynaLearn(); % ~ 1sec
    m = m.dlLoad(char("models/dlPredictiveCorticalCircuitModelLWK" + string(i))); % ~ 10sec, New larger model; keeping track of its activity in Gamma/Beta **
   
    ratioLog1 = {m.dlCustomLog(1, :)};
    ratioLog1 = cell2mat(ratioLog1{1});
    ratioLogs1{i} = ratioLog1';

    ratioLog2 = {m.dlCustomLog(2, :)};
    ratioLog2 = cell2mat(ratioLog2{1});
    ratioLogs2{i} = ratioLog2';

    ratioLog3 = {m.dlCustomLog(3, :)};
    ratioLog3 = cell2mat(ratioLog3{1});
    ratioLogs3{i} = ratioLog3';

    ratioLog4 = {m.dlCustomLog(4, :)};
    ratioLog4 = cell2mat(ratioLog4{1});
    ratioLogs4{i} = ratioLog4';

    lfpcl = {m.dlCustomLog(5, :)};
    lfps(i, :, :) = cell2mat(lfpcl{1});

end

%% To matrix

clc;

ratioLogs1 = cell2mat(ratioLogs1);
ratioLogs2 = cell2mat(ratioLogs2);
ratioLogs3 = cell2mat(ratioLogs3);
ratioLogs4 = cell2mat(ratioLogs4);

ratioLogs1 = ratioLogs1'; % Theta
ratioLogs2 = ratioLogs2'; % Alpha
ratioLogs3 = ratioLogs3'; % Beta
ratioLogs4 = ratioLogs4'; % Gamma

%% Averages

clc;
ratioLogsPooled = ratioLogs1 / mean(ratioLogs1, 'all');
ratioLogsAvg = zeros(size(ratioLog1));

for i = 1:150

    ratioLogsPooled(i, :) = conv(ratioLogsPooled(i, :), [1 1 1]/3, "same");
%     ratioLogsPooled(i, :) = ratioLogsPooled(i, :) - mean(ratioLogsPooled(i, :));

end

for i = 1:15

    q = ratioLogsPooled(i:10:150, :);
    ratioLogsAvg(i, :) = mean(q, 1);
%     ratioLogsAvg(i, :) = ratioLogsAvg(i, :) - ratioLogsAvg(i, 5);

end

%% Ratio change (relative) in trials during predictive task training

clc;

w = 100;
cnt = 1;
tlab = ["A", "U1", "B", "U2", "C", "U3", "A", "U4", "B", "U5", "C", "U6"];
figure("Position", [0 0 1500 1000]);

set1 = 2:9;
% ratioLogPlot = ratioLogsAvg(:, 5:590);

for i = 1:8

    for l = 0:5

        subplot(4, 2, i);
        fill([l*w, l*w, l*w+w, l*w+w], [-1e+5, 1e+5, 1e+5, -1e+5], [sin(l*0.2), 1, cos(l*0.2)], 'HandleVisibility','off');
        hold('on');

    end

end

dispLabels = m.dlCustomLogLabel{1};

for i = set1

    subplot(4, 2, mod(cnt-1, 8)+1);
    legend("Location", "southwest");

    x = 5:590;
    q = ratioLogsPooled(i:10:150, x);
                                          % Create Independent Variable
    y = q;                                  % Create Dependent Variable 6Experimentsz Data
    N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
    yMean = mean(y, 1);                                    % Mean Of All Experiments At Each Value Of ,x 
    ySEM = 2.14*std(y)/sqrt(N);                              % Compute 2Standard Error Of The Mean* Of All Experiments At Each Value Of .x1
    y1 = yMean + ySEM;
    y2 = yMean - ySEM;
%     CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
%     yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

    ql = min([y1, y2]);
    qu = max([y1, y2]);

    for l = 0:5
        text((l+0.5)*w, qu+.1, tlab(l+1));
    end

    ylim([ql-.2 qu+.2]);
    title(dispLabels{i});
    
    plot(x, yMean, '-', 'DisplayName', dispLabels{i});                                 % Plot Mean Of All Experiments
    hold("on");
%     plot(x, yCI95+yMean,'-r')                           % Plot 95% Confidence Intervals Of All Experiments
    patch([x, fliplr(x)], [y1, fliplr(y2)], 'b', 'EdgeColor','g', 'FaceAlpha',0.25, 'DisplayName', 'CI-95%');

    hold("on");
    grid("on");
    cnt = cnt + 1;

end

sgtitle("Theta power band ratio change (scaled and relative) across trials (Average across 10 different models)");

%% Ratio change (relative) in trials during predictive task training

clc;

w = 100;
cnt = 1;
tlab = ["P", "U", "P", "U", "P", "U"];
f = figure("Position", [0 0 1500 1000]);

set1 = 2:9;

ratioLogPlot = ratioLogsAvg(:, 5:590);

for i = 1:8

    for l = 2:3

        subplot(4, 2, i);
        fill([l*w, l*w, l*w+w, l*w+w], [-40, 40, 40, -40], [sin(l*0.35), 1, cos(l*0.35)], 'HandleVisibility','off');
        hold('on');

    end

end

dispLabels = m.dlCustomLogLabel{1};

for i = set1

    subplot(4, 2, mod(cnt-1, 8)+1);
%     plot(ratioLogPlot(i, :), 'DisplayName', dispLabels{i});
    legend("Location", "southwest");

%     for l = 0:5
%         text((l+0.5)*w, qu+0.5, tlab(l+1));
%     end

    x = 200:400;
    q = ratioLogsPooled(i:10:150, x);
                                          % Create Independent Variable
    y = q;                                  % Create Dependent Variable 6Experimentsz Data
    N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
    yMean = mean(y, 1);                                    % Mean Of All Experiments At Each Value Of ,x 
    ySEM = 2*std(y)/sqrt(N);                              % Compute 2Standard Error Of The Mean* Of All Experiments At Each Value Of ‘x’
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

sgtitle("Theta power band ratio change from (Average across 10 different models/ Predictable block (green) to unpredictable (yellow))");

%%

clc;
for i = 1

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
trial = 20;
k = 5;

lfptemp = {m.dlCustomLog{5, :}};
lfptemp = cell2mat(lfptemp);

%%

figure("Position", [0 0 1500 1000]);
subplot(1, 1, 1);

for i = 2:popsize

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

%% End
