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


% obj = obj.squeeze;
pop_axis = obj.findaxis('populations'); if isempty(pop_axis); error('obj must contain a populations axis to be converted to DynaSim'); end
var_axis = obj.findaxis('variables'); if isempty(var_axis); error('obj must contain a variables axis to be converted to DynaSim'); end
varied_inds = true(1,ndims(obj)); varied_inds(pop_axis) = false; varied_inds(var_axis) = false;
varied_axis = find(varied_inds);

% Bring pop and var to front
obj = obj.permute([pop_axis,var_axis, varied_axis(:)']);

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
    
    
end

data = CheckData(data);

end


function varied_names = only_varieds(all_names)
    inds = true(1,length(all_names));
    inds(strcmp(all_names,'populations')) = false; 
    inds(strcmp(all_names,'variables')) = false;
    varied_names = all_names(inds);
end

