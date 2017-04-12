function data_out = downsampleData(data,ds)
%DOWNSAMPLEDATA - Downsamples DynaSim data structre data
%
% Usage:
%   data_out = ds.calcFR(data,'option',value)
%
% Inputs:
%   - data: DynaSim data structure (see ds.checkData)
%   - ds: number of datapoints to downsample
%
% Outputs:
%   - data_out: data structure all fields replaced by their average valeus
%               (averaged across neurons).
%
% See also: ds.plotFR, ds.analyzeStudy, dsSimulate, ds.checkData, ds.selectVariables

%% 1.0 Check inputs

data = ds.checkData(data, varargin{:});
% note: calling ds.checkData() at beginning enables analysis function to
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
