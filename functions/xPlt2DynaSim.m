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


pop_axis = obj.findaxis('populations')
obj = obj.mergeDims([1:pop_axis-1]);
obj = obj.mergeDims([pop_axis:ndims(obj)]);
obj = squeeze(obj); % Squeeze out the empty dimensions.
% obj.getaxisinfo;

data = struct;
ax_vals = obj.exportAxisVals;
ax_names = obj.exportAxisNames;
varied = obj.axis(1).astruct.premerged_names;
varied_vals = obj.axis(1).astruct.premerged_values;
for i = 1:size(obj,1)
    % Add actual data
    for j = 1:size(obj,2)
        data(i).(ax_vals{2}{j}) = obj.data{i,j};
    end
    
    % Add list of varied variables
    data(i).varied = varied;
    
    % Add values of varied variables
    for j = 1:length(varied)
        data(i).(varied{j}) = varied_vals{j}(i);
    end
    
    % Add other DynaSim info if present
    obj_curr = obj.meta.dynasim;
    fc = 'labels'; if isfield(obj_curr,fc); data(i).(fc) = obj_curr.(fc); end
    fc = 'model'; if isfield(obj_curr,fc); data(i).(fc) = obj_curr.(fc); end
    fc = 'simulator_options'; if isfield(obj_curr,fc); data(i).(fc) = obj_curr.(fc); end
    fc = 'time'; if isfield(obj_curr,fc); data(i).(fc) = obj_curr.(fc); end
    
    
end


end
