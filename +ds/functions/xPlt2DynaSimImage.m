function data=xPlt2DynaSimImage(obj)
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
% d2b = xPlt2DynaSim(squeeze(DynaSim2xPlt(data2)));
% % Make sure 1 is identical
% close all; 
% PlotData(data1); PlotData(d1);
% % Make sure 2 is identical
% PlotData(data2); PlotData(d2);
% PlotData(data2); PlotData(d2b);

% Squeeze out populations and variables axes (should remove them if they
% are of size 1)
obj = obj.squeezeRegexp('populations');
obj = obj.squeezeRegexp('variables');

% Find population and variable axes
if ~isempty(obj.findaxis('populations')); error('Populations axis should be empty'); end
if ~isempty(obj.findaxis('variables')); error('Variables axis should be empty'); end

% Merge all varieds together
obj = obj.mergeDims(1:ndims(obj));

% Get rid of any leftover axes created by mergeDims
obj = obj.squeezeRegexp('Dim');

% Make sure the "varied" parameters are along dimension 2 (e.g. so it's
% 1xNvaried)
if iscolumn(obj.data); obj=obj.transpose; end

% Build DynaSim data structure
data = struct;
ax_vals = obj.exportAxisVals;
ax_names = obj.exportAxisNames;
varied = obj.axis(2).astruct.premerged_names;
varied_vals = obj.axis(2).astruct.premerged_values;
has_varied = 1;

for j = 1:size(obj,2)                               % Loop through varieds
    
    data(j).plot_files = obj.data{j};
    
    % If there are any varied parameters....
    if has_varied
        % Add list of varied variables
        data(j).varied = varied;

        % Add values of varied variables
        for i = 1:length(varied)
            data(j).(varied{i}) = varied_vals{i}(j);
        end
    end
    
end


end
