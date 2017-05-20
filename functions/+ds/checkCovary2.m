function [effective_vary_lengths, linked_indices] = checkCovary2(vary_lengths, vary_params, data_length)
  %CHECKCOVARY - TODO I assume this checks if any varied parameters are covaried?

  % Remove any vary statements used to set parameter values.
  non_singleton_vary_lengths = vary_lengths(vary_lengths > 1);
  n_varied = length(non_singleton_vary_lengths);

  %% Check all combinations of varied parameters to see if the product of the number
  % of their unique values describes the amount of data simulated.
  full_rank_choices = {};
  index = 1;

  for n_chosen = 1:n_varied
    % Check all pairs, then all triplets...
    choices = nchoosek(1:n_varied, n_chosen);
    n_choices = nchoosek(n_varied, n_chosen);

    % For each subset of varied parameters, see if the product of the numbers of
    % unique parameters equals data_length.
    for choice = 1:n_choices
      estimated_data_length(index) = prod(non_singleton_vary_lengths(choices(choice, :)));
      if estimated_data_length == data_length;
        full_rank_choices{index} = choices(choice, :);
        index = index + 1;
      end
    end

  end

  % Remove linked sets of parameters from effective_vary_lengths.
  marked_for_removal = zeros(size(effective_vary_lengths));
  for l = 1:number_linked_sets
    marked_for_removal(linked_indices{l}(2:end)) = 1;
  end
  effective_vary_lengths = non_singleton_vary_lengths;
  effective_vary_lengths(marked_for_removal) = [];

end

function linked_indices = find_linked_params(full_rank_choices)
  %% Parse the results of finding "full rank" choices.
  if isempty(full_rank_choices)
    % No combination of varied parameters equals total data length,
    % so treat data as linear.
    linked_indices{1} = 1:length(non_singleton_vary_lengths);

  elseif length(full_rank_choices) == 1
    % Only one combination of varied parameters describes the total data length,
    % so there are no covaried parameters.
    linked_indices = {};

  else
    % If more than one combination of varied parameters describes the total data length,
    % there are covaried parameters, and we have to find them.
    no_full_rank_choices = length(full_rank_choices);

    % Create a matrix in which each row represents a full rank choice of varied parameters.
    full_rank_participation = zeros(no_full_rank_choices, n_varied);
    for frc = 1:no_full_rank_choices
      full_rank_participation(frc, full_rank_choices{frc}) = 1;
    end

    % If a certain parameter is not involved in any "full rank" product, remove
    % it from effective_vary_lengths.
    non_participating_parameters = sum(full_rank_participation) == 0;
    effective_vary_lengths(non_participating_parameters) = [];

    % For each pair of full rank choices, see if they differ in exactly two parameters.
    frc_pairs = nchoosek(1:no_full_rank_choices, 2);
    no_frc_pairs = size(frc_pairs, 1);
    index = 1;
    for pair = 1:no_frc_pairs
      frc_diff = full_rank_choices(frc_pairs(pair, 1), :) - full_rank_choices(frc_pairs(pair, 2), :);
      if sum(abs(frc_diff) > 0) == 2
        linked_indices{index} = find(abs(frc_diff) == 1);
        index = index + 1;
      else
        warning(' dimensions of parameter set ambiguous; data will be plotted as linear.')
        clear linked_indices
        linked_indices{1} = 1:length(non_singleton_vary_lengths);
        return
    end

  end

end
