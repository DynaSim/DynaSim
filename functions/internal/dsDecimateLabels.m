function data = dsDecimateLabels(data,labels2keep, varargin)
% dsDecimateData - Downsamples DynaSim data structure by removing all but a
% subset of labels
%
% Usage:
%   data_out = dsDecimateData(data,labels2keep,'option',value)
%
% Inputs:
%   - data: DynaSim data structure (see dsCheckData)
%   - labels2keep: a subset of data.labels that describes the labels to
%   keep
%
% Outputs:
%   - data: data structure all fields replaced by their decimated values
%
% See also: dsDecimateCells, dsDecimateData

%% 1.0 Check inputs

data = dsCheckData(data, varargin{:});
% note: calling dsCheckData() at beginning enables analysis function to
% accept data matrix [time x cells] in addition to DynaSim data structure.

%% do the decimating    
% Loop through all labels and remove all but the ones that are in labels2keep
for j = 1:length(data(1).labels)

    % Is this a label to keep?
    tokeep = any(strcmp(data(1).labels{j}, labels2keep));

    % If not, remove from structure
    if ~tokeep
        data = rmfield(data,data(1).labels{j});
    end

end



%%
% More condensed version
data_decim = data;

allfields = fieldnames(data_decim);

data_decim2 = [];
for i = 1:length(data_decim)
    dpart = struct;
    for j = 1:length(allfields)
        dpart.(allfields{j}) = data_decim(i).(allfields{j});
    end
    data_decim2 = [data_decim2, dpart];
end

