function [effective_vary_lengths, linked_indices] = CheckCovary2(vary_lengths, vary_params, data_length)
%CHECKCOVARY - TODO I assume this checks if any varied parameters are covaried?

non_singleton_vary_lengths = vary_lengths(vary_lengths > 1);

n_varied = length(vary_lengths);

full_rank_choices = {}; % full_rank = nan(n_varied*2^(n_varied - 1), 1);

index = 1;

for n_chosen = 1:n_varied
    
    choices = nchoosek(1:n_varied, n_chosen);
    
    n_choices = nchoosek(n_varied, n_chosen);
   
    for choice = 1:n_choices
        
        estimated_data_length(index) = prod(non_singleton_vary_lengths(choices(choice, :)));
        
        index = index + 1;
        
    end
    
end

if isempty(full_rank_choices)
    
    % Check all pairs, then all triplets... Or, to optimize, check sum of
    % all, then each and sums of n - 1, then pairs and sums of n - 2,
    % etc....
    
    linked_indices{1} = 1:length(non_singleton_vary_lengths);
    
end
    
    effective_vary_lengths = non_singleton_vary_lengths;
    
    for l = 1:number_linked_sets
        
        effective_vary_lengths(linked_indices{l}(2:end)) = [];
        
    end

    effective_vary_lengths = non

full_rank_choices{:}