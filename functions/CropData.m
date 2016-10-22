function data_out = CropData(data,inds)
%% data_out = CropData(data,'option',value)
% Purpose: Crops DynaSim data structre
% Inputs:
%   data - DynaSim data structure (see CheckData)
%   inds - data points to retain
% Outputs:
%   data_out: data structure with all data cropped 


%% 1.0 Check inputs

data = CheckData(data);
% note: calling CheckData() at beginning enables analysis function to
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