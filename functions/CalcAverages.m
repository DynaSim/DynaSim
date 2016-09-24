function data_out = CalcAverages(data)
%% data_out = CalcFR(data,'option',value)
% Inputs:
%   data - DynaSim data structure (see CheckData)
% Outputs:
%   data_out: data structure all fields replaced by their average valeus
%             (averaged across neurons).

% See also: PlotFR, AnalyzeStudy, SimulateModel, CheckData, SelectVariables

%% 1.0 Check inputs

data = CheckData(data);
% note: calling CheckData() at beginning enables analysis function to
% accept data matrix [time x cells] in addition to DynaSim data structure.

%% do the averaging

data_out = data;
for i = 1:length(data)
    % Identify all fields in data containing simulated output
    labels = data(i).labels;
    labels = labels(~strcmp(labels,'time'));
    
    % Sweep through these fields and take average
    for j = 1:length(labels)
        data_out(i).(labels{j}) = mean(data(i).(labels{j}),2);
    end
    
    % Set the size of all populations to 1
    for j = 1:length(data(i).model.specification.populations)
        data_out(i).model.specification.populations(j).size = 1;
    end
    
end






end