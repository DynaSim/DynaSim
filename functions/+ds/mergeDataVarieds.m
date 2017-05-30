
function [data_new, variedname_merged, varied_vals ] = mergeDataVarieds(data,varied_fields)
    % [data_new, variedname_merged, varied_vals ] = mergeDataVarieds(data,varied_fields)
    %
    % Purpose: This function takes in the current DynaSim datastructure, data, and
    % returns data_new. data_new is the same as the original, except the
    % fields listed in varied_fields are merged together. This is useful
    % for 
    %
    % Usage:
    %   [data_new, variedname_merged, varied_vals ] = mergeDataVarieds(data,varied_fields)
    %
    % Inputs:
    %   data: DynaSim data structure
    %   varied_fields: cell array of field names that will be merged
    %   together
    %
    % Outputs:
    %   data_new: DynaSim output structure
    %   variedname_merged: New varied fieldname assigned to the merge varied fields
    %   varied_vals: New values assigned to the merge varied fields
    %
    % Examples:
    %   [data, variedname_merged, varied_vals ] = mergeVarieds(data,vary_labels(linked_inds{j}));
    %   (see ds2mdd for example)
    % 
    % Submodules: cat_with_underscores
    %
    % Author: David Stanley, Boston University, 2017
    %
    % See also: unmergeDataVarieds, modifications2Vary, vary2Modifications
    % 

    varied_vals = cell(length(data),length(varied_fields));
    for i = 1:length(data)
        data_temp = data(i);
        variedname_merged = cat_with_underscores(varied_fields);
        variedname_merged = strcat('Covaried_',variedname_merged);
        
        for j = 1:length(varied_fields)
            
            varied_vals{i,j} = data_temp.(varied_fields{j});       % Store varied value. Storead as varied x sim number
            data_temp = rmfield(data_temp,varied_fields{j});         % Remove original fields
        end
        
        varied_str{i} = cellfun(@(x) num2str(x),varied_vals(i,:),'UniformOutput',0);
        varied_str{i} = cat_with_underscores(varied_str{i});
        data_temp.(variedname_merged) = varied_str{i};
        
        [~, inds] = intersect(data_temp.varied,varied_fields);
        inds2 = true(1,length(data_temp.varied));
        inds2(inds) = false;
        data_temp.varied = data_temp.varied(inds2);
        data_temp.varied{end+1} = variedname_merged;
        
        data_new(i) = data_temp;
    end
    
    
    

end

function str_out = cat_with_underscores(cellstr_in)
% Takes in a cell array of chars and concatenates them together with
% underscores separating the original divisions between cells. E.g.
% {'cat','dog'} becomes 'cat_dog'
    
    temp = vertcat(cellstr_in(:)', repmat({'_'},1,length(cellstr_in)));
    temp = temp(:)';
    str_out = horzcat(temp{1:end-1});
end