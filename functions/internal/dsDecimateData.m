function data_out = dsDecimateData(data,ds, varargin)
%DOWNSAMPLEDATA - Downsamples DynaSim data structre data
%
% Usage:
%   data_out = dsCalcFR(data,'option',value)
%
% Inputs:
%   - data: DynaSim data structure (see dsCheckData)
%   - ds: number of datapoints to downsample
%
% Outputs:
%   - data_out: data structure all fields replaced by their decimated values
%
% See also: dsPlotFR, dsAnalyzeStudy, dsSimulate, dsCheckData, dsSelectVariables

%% 1.0 Check inputs

data = dsCheckData(data, varargin{:});
% note: calling dsCheckData() at beginning enables analysis function to
% accept data matrix [time x cells] in addition to DynaSim data structure.

%% do the decimating
data_out = data;
for i = 1:length(data)
    % Identify all fields in data containing simulated output
    labels = data(i).labels;
    
    % Sweep through these fields and decimate
    for j = 1:length(labels)
        data_out(i).(labels{j}) = data(i).(labels{j})(1:ds:end,:);
    end
    
end


end
