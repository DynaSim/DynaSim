function data_out = cropData(data,inds)
%CROPDATA - Crops DynaSim data structre
%
% Usage:
%   data_out = ds.cropData(data,'option',value)
%
% Inputs:
%   - data: DynaSim data structure (see ds.checkData)
%   - inds: data points to retain
%
% Outputs:
%   - data_out: data structure with all data cropped

%% 1.0 Check inputs

data = ds.checkData(data);
% note: calling ds.checkData() at beginning enables analysis function to
% accept data matrix [time x cells] in addition to DynaSim data structure.

%% do the cropping

data_out = data;
for i = 1:length(data)
    % Identify all fields in data containing simulated output
    labels = data(i).labels;

    % Sweep through these fields and take average
    for j = 1:length(labels)
        data_out(i).(labels{j}) = data(i).(labels{j})(inds,:);
    end
end

end
