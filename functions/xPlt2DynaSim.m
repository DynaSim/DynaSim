function data=xPlt2DynaSim(obj)
%% data=xPlt2DynaSim(obj,varargin)
% Dependencies:
%   Requires the MDD class, which should be part of DynaSim. If not,
%   get it here https://github.com/davestanley/MDD
% Author: David Stanley, Boston University, stanleyd@bu.edu

% This combines the 2D "vary" sweep into a single dimension. It also
% combines all populations and variables into a single 1D list. Thus, Axis
% 1 is equivalent to Jason's structure array - data(1:9). Axis 4 is
% equivalent to the structure fields in Jason's DynaSim structure.
% 
% % % Testing code
% load sample_data_dynasim.mat
% data1=data;
% data2 = data(1);
% d1 = xPlt2DynaSim(DynaSim2xPlt(data1));
% d2 = xPlt2DynaSim(DynaSim2xPlt(data2));
% % Make sure 1 is identical
% close all; PlotData(data1); PlotData(d1);
% % Make sure 2 is identical
% close all; PlotData(data2); PlotData(d2);


% obj = obj.squeeze;
pop_axis = obj.findaxis('populations');
var_axis = obj.findaxis('variables');
varied_inds = true(1,ndims(obj)); varied_inds(pop_axis) = false; varied_inds(var_axis) = false;
varied_axis = find(varied_inds);

% Bring pop and var to front
obj = obj.permute([pop_axis,var_axis, varied_axis(:)']);

% Find population sizes for each population (will be used MUCH later)
num_pops = size(obj,1);
pop_names = obj.exportAxisVals; pop_names = pop_names{1};
population_sz_all = cellfun(@(x) size(x,2),obj.data,'UniformOutput',1);


% Merge populations and variables together
obj = obj.mergeDims(1:2);

% Merge all other varied variables
obj = obj.mergeDims(3:ndims(obj));

% Build DynaSim data structure
data = struct;
ax_vals = obj.exportAxisVals;
ax_names = obj.exportAxisNames;
varied = obj.axis(3).astruct.premerged_names;
varied_vals = obj.axis(3).astruct.premerged_values;

obj = obj.squeeze;
for j = 1:size(obj,2)
    
    % Add actual data
    for i = 1:size(obj,1)
        data(j).(ax_vals{1}{i}) = obj.data{i,j};
    end
    
    % If there are any varied parameters....
    if 1
        % Add list of varied variables
        data(j).varied = varied;

        % Add values of varied variables
        for i = 1:length(varied)
            data(j).(varied{i}) = varied_vals{i}(j);
        end
    end
    
    % Add other DynaSim info if present
    obj_curr = obj.meta.dynasim;
    fc = 'labels'; if isfield(obj_curr,fc); data(j).(fc) = obj_curr.(fc); end
    fc = 'model'; if isfield(obj_curr,fc); data(j).(fc) = obj_curr.(fc); end
    fc = 'simulator_options'; if isfield(obj_curr,fc); data(j).(fc) = obj_curr.(fc); end
    fc = 'time'; if isfield(obj_curr,fc); data(j).(fc) = obj_curr.(fc); end
    
    
    for i = 1:length(data)
        for j = 1:num_pops
            obj_temp = obj(['^' pop_names{j}],i);   % ^ is regexp for begins with
            
            % Get size of all variables in this population
            pop_sz = cellfun(@(x) size(x,2),obj_temp.data);
            pop_sz2 = pop_sz(pop_sz ~= 0);  % Get rid of empties if any
            
            % Make sure that all variables have the same population size
            if any(pop_sz2 ~= mode(pop_sz2)); error('Variables in population %s have differing population sizes',pop_names{j}); end
            
            % Assign population size to model info
            data(i).model.specification.populations(j).size = mode(pop_sz2);
        end
        data(i).labels = {ax_vals{1}{:}, 'time'};
    end
    
    
end

data = CheckData(data);

end


function varied_names = only_varieds(all_names)
    inds = true(1,length(all_names));
    inds(strcmp(all_names,'populations')) = false; 
    inds(strcmp(all_names,'variables')) = false;
    varied_names = all_names(inds);
end

