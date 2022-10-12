%% Load model logs (B/G ratio)

clear;clc;

ratioLogs = cell(10, 1);

for i = 1:10

    m = DynaLearn(); % ~ 1sec
    m = m.dlLoad(char("models/dlPredictiveCorticalCircuitModelLWK" + string(i))); % ~ 10sec, New larger model; keeping track of its activity in Gamma/Beta **
   
    ratioLog = {m.dlCustomLog(1, :)};
    ratioLog = cell2mat(ratioLog{1});
    ratioLogs{i} = ratioLog;

end

%% To matrix

clc;
ratioLogs = cell2mat(ratioLogs);

%% Averages

ratioLogsPooled = ratioLogs;
ratioLogsAvg = zeros(size(ratioLog));

for i = 1:150

    ratioLogsPooled(i, :) = conv(ratioLogs(i, :), [1 1 1]/3, "same");
    ratioLogsPooled(i, :) = ratioLogsPooled(i, :) - mean(ratioLogsPooled(i, :));

end

for i = 1:15

    q = ratioLogsPooled(i:10:150, :);
    ratioLogsAvg(i, :) = mean(q, 1);
    ratioLogsAvg(i, :) = ratioLogsAvg(i, :) - ratioLogsAvg(i, 5);

end

%% Beta to Gamma ratio change (relative) in trials during predictive task training

clc;

w = 50;
cnt = 1;
tlab = ["A", "U1", "B", "U2", "C", "U3", "A", "U4", "B", "U5", "C", "U6"];
f = figure("Position", [0 0 1500 1000]);

set1 = 2:9;

% ratioLog = {m.dlCustomLog(1, :)};
% ratioLog = cell2mat(ratioLog{1});
ratioLogPlot = ratioLogsAvg(:, 5:290);

for i = 1:8

    for l = 0:5

        subplot(4, 2, i);
        fill([l*w, l*w, l*w+w, l*w+w], [-40, 40, 40, -40], [sin(l*0.2), 1, cos(l*0.2)], 'HandleVisibility','off');
        hold('on');

    end

end

dispLabels = m.dlCustomLogLabel{1};

for i = set1

    subplot(4, 2, mod(cnt-1, 8)+1);
%     plot(ratioLogPlot(i, :), 'DisplayName', dispLabels{i});
    legend("Location", "southwest");

    for l = 0:5
        text((l+0.5)*w, qu+0.5, tlab(l+1));
    end

    x = 5:290;
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
    title(dispLabels{i});
    
    plot(x, yMean, '-')                                  % Plot Mean Of All Experiments
    hold("on");
%     plot(x, yCI95+yMean,'-r')                           % Plot 95% Confidence Intervals Of All Experiments
    patch([x, fliplr(x)], [y1, fliplr(y2)], 'b', 'EdgeColor','g', 'FaceAlpha',0.25)

    hold("on");
    grid("on");
    cnt = cnt + 1;

end

sgtitle("Relative Beta to Gamma power band ratio change from first trial (Average across 10 different models)");

%% Beta to Gamma ratio change (relative) in trials during predictive task training

clc;

w = 50;
cnt = 1;
tlab = ["P", "U", "P", "U", "P", "U"];
f = figure("Position", [0 0 1500 1000]);

set1 = 2:9;

ratioLogPlot = ratioLogsAvg(:, 5:290);

for i = 1:8

    for l = 0:5

        subplot(4, 2, i);
        fill([l*w, l*w, l*w+w, l*w+w], [-40, 40, 40, -40], [sin(l*0.3), 1, cos(l*0.3)], 'HandleVisibility','off');
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

    x = 100:200;
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
    
    plot(x, yMean, '-')                                  % Plot Mean Of All Experiments
    hold("on");
%     plot(x, yCI95+yMean,'-r')                           % Plot 95% Confidence Intervals Of All Experiments
    patch([x, fliplr(x)], [y1, fliplr(y2)], 'b', 'EdgeColor','g', 'FaceAlpha',0.25)

    hold("on");
    grid("on");
    cnt = cnt + 1;

end

sgtitle("Relative Beta to Gamma power band ratio change from first trial (Average across 10 different models/ Predictable block (blue) to unpredictable (yellow))");


%% End