function data_out = CalcAverages(data,opts)
%% data_out = CalcAverages(data,'option',value)
% Inputs:
%   data - DynaSim data structure (see CheckData)
%   opts - Options structure
%          opts.save_std - Also creates new fields storing standard deviations (boolean; true/false)

% Outputs:
%   data_out: data structure all fields replaced by their average vales
%             (averaged across neurons).

% See also: PlotFR, AnalyzeStudy, SimulateModel, CheckData, SelectVariables

%% 1.0 Check inputs

if nargin < 2
    opts = struct;
end

data = CheckData(data);
% note: calling CheckData() at beginning enables analysis function to
% accept data matrix [time x cells] in addition to DynaSim data structure.

opts = struct_addDef(opts,'save_std',0);

%% do the averaging

data_out = data;
for i = 1:length(data)
    % Identify all fields in data containing simulated output
    labels = data(i).labels;
    labels = labels(~strcmp(labels,'time'));
    
    % Sweep through these fields and take averages
    for j = 1:length(labels)
        % Save mean
        data_out(i).(labels{j}) = mean(data(i).(labels{j}),2);
        
        % Save standard deviation also if requested
        if opts.save_std
            data_out(i).([labels{j} '_std']) = std(data(i).(labels{j}),0,2);
        end
    end
    
    % Set the size of all populations to 1
    for j = 1:length(data(i).model.specification.populations)
        data_out(i).model.specification.populations(j).size = 1;
    end
    
end



end