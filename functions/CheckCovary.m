function [effective_vary_lengths, linked_indices] = CheckCovary(vary_lengths, vary_params)

non_singleton_vary_lengths = vary_lengths(vary_lengths > 1);

diff_vary_lengths = diff(non_singleton_vary_lengths);

equal_vary_length_indices = diff_vary_lengths == 0;

equal_vary_length_blocks = IndexToBlocks(equal_vary_length_indices);

if ~isempty(equal_vary_length_blocks)

equal_vary_length_blocks(:, 2) = equal_vary_length_blocks(:, 2) + 1;

end

no_blocks = size(equal_vary_length_blocks, 1);

linked_indices = {};
    
number_linked_sets = 0;

for block = 1:no_blocks
    
    param_indices = equal_vary_length_blocks(block, 1):equal_vary_length_blocks(block, 2);
    
    iteration_number = 1;
    
    while length(param_indices) > 1
        
        [linked_set, param_indices] = find_linked_params(param_indices, vary_params);
        
        if ~isempty(linked_set)
            
            number_linked_sets = number_linked_sets + 1;
            
            linked_indices{number_linked_sets} = linked_set;
            
        end
        
        iteration_number = iteration_number + 1;
        
        if iteration_number > 10
            
            return
            
        end
        
    end
    
end

effective_vary_lengths = vary_lengths;

for l = 1:number_linked_sets
    
    effective_vary_lengths(linked_indices{l}(2:end)) = [];
    
end

end

function [linked_indices, non_linked_indices] = find_linked_params(param_indices, vary_params)

param1_index = param_indices(1);

params_at_value1 = vary_params(vary_params(:, param1_index) == vary_params(1, param1_index), :);

no_values_at_value1 = nan(size(param_indices));

no_values_at_value1(1) = 1;

for p = 2:length(param_indices)
    
    param_index = param_indices(p);
    
    no_values_at_value1(p) = length(unique(params_at_value1(:, param_index)));
    
end

linked_indices = param_indices(no_values_at_value1 == 1);

non_linked_indices = param_indices(no_values_at_value1 > 1);

end