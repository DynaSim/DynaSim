function data_out = DownsampleData(data,ds)
%DOWNSAMPLEDATA - Downsamples DynaSim data structre data
%
% Usage:
%   data_out = calcFR(data,'option',value)
%
% Inputs:
%   - data: DynaSim data structure (see checkData)
%   - ds: number of datapoints to downsample
%
% Outputs:
%   - data_out: data structure all fields replaced by their average valeus
%               (averaged across neurons).
%
% See also: PlotFR, AnalyzeStudy, SimulateModel, checkData, SelectVariables

%% 1.0 Check inputs

data = checkData(data);
% note: calling checkData() at beginning enables analysis function to
% accept data matrix [time x cells] in addition to DynaSim data structure.

%% do the averaging

data_out = data;
for i = 1:length(data)
    % Identify all fields in data containing simulated output
    labels = data(i).labels;
    
    % Sweep through these fields and take average
    for j = 1:length(labels)
        data_out(i).(labels{j}) = data(i).(labels{j})(1:ds:end,:);
    end
    
end






end
