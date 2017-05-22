function [effective_vary_indices, linked_indices] = checkCovary2(vary_lengths, data_length, varargin)
  %CHECKCOVARY - TODO I assume this checks if any varied parameters are covaried?
 
  %% localfn output
  if ~nargin
      effective_vary_lengths = localfunctions; % output var name specific to this fn
      return
  end
  
  %% auto_gen_test_data_flag argin
  options = ds.checkOptions(varargin,{'auto_gen_test_data_flag',0,{0,1}},false);
  if options.auto_gen_test_data_flag
      varargs = varargin;
      varargs{find(strcmp(varargs, 'auto_gen_test_data_flag'))+1} = 0;
      varargs(end+1:end+2) = {'unit_test_flag',1};
      argin = [{vary_lengths}, {data_length}, varargs]; % specific to this function
  end

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
      estimated_data_length = prod(non_singleton_vary_lengths(choices(choice, :)));
      if estimated_data_length == data_length
        full_rank_choices{index} = choices(choice, :);
        index = index + 1;
      end
    end

  end
  
  linked_indices = find_linked_params(full_rank_choices, length(non_singleton_vary_lengths), varargin{:});

  % Remove linked sets of parameters from effective_vary_lengths.
  marked_for_removal = zeros(size(non_singleton_vary_lengths));
  for l = 1:length(linked_indices)
    marked_for_removal(linked_indices{l}(2:end)) = 1;
  end
  % effective_vary_lengths = non_singleton_vary_lengths;
  % effective_vary_lengths(logical(marked_for_removal)) = [];
  effective_vary_indices = ~logical(marked_for_removal);
  
  %% auto_gen_test_data_flag argout
  if options.auto_gen_test_data_flag
      argout = {effective_vary_lengths, linked_indices}; % specific to this function
      
      ds.unit.saveAutoGenTestData(argin, argout);
  end

end

function linked_indices = find_linked_params(full_rank_choices, no_varied, varargin)

  %% auto_gen_test_data_flag argin
  options = ds.checkOptions(varargin,{'auto_gen_test_data_flag',0,{0,1}},false);
  if options.auto_gen_test_data_flag
    varargs = varargin;
    varargs{find(strcmp(varargs, 'auto_gen_test_data_flag'))+1} = 0;
    varargs(end+1:end+2) = {'unit_test_flag',1};
    argin = [{full_rank_choices}, varargs];
  end

  linked_indices = {};
  
  %% Parse the results of finding "full rank" choices.
  if isempty(full_rank_choices)
    % No combination of varied parameters equals total data length,
    % so treat data as linear.
    linked_indices{1} = 1:no_varied;

  else
    % If one or more combinations of varied parameters describes the total data length,
    % we need to link covaried parameters.
    no_full_rank_choices = length(full_rank_choices);

    % Create a matrix in which each row represents a full rank choice of varied parameters.
    full_rank_participation = zeros(no_full_rank_choices, no_varied);
    for frc = 1:no_full_rank_choices
      full_rank_participation(frc, full_rank_choices{frc}) = 1;
    end
    
    % If a certain parameter is not involved in any "full rank" product,
    % consider it linked.
    index = 1;
    non_participating_parameters = find(sum(full_rank_participation, 1) == 0);
    if ~isempty(non_participating_parameters)
        linked_indices{index} = [0 non_participating_parameters];
        index = index + 1;
    end

    if no_full_rank_choices > 1
        % If more than one combination of varied parameters describes the total
        % data length, for each pair of full rank choices, see if they differ 
        % in exactly two parameters.
        frc_pairs = nchoosek(1:no_full_rank_choices, 2);
        no_frc_pairs = size(frc_pairs, 1);
        for pair = 1:no_frc_pairs
            frc_diff = full_rank_participation(frc_pairs(pair, 1), :) - full_rank_participation(frc_pairs(pair, 2), :);
            if sum(abs(frc_diff) > 0) == 2
                linked_indices{index} = find(abs(frc_diff) == 1);
                index = index + 1;
                %   else
                %     warning(' dimensions of parameter set ambiguous; data will be plotted as linear.')
                %     clear linked_indices
                %     linked_indices{1} = 1:no_varied;
                %     return
            end
        end
        
    end

  end
  
  %% auto_gen_test_data_flag argout
  if options.auto_gen_test_data_flag
      argout = {linked_indices}; % specific to this function
      
      ds.unit.saveAutoGenTestDataLocalFn(argin, argout); % localfn
  end
  
end
