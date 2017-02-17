

function [ind_dim,ind_embedded,ind_selected] = xPltSplice (s,dim_name,selection_criteria,string_matching_func)

    if nargin < 4
        string_matching_func = @regexp;
    end
    
    if nargin < 3
        selection_criteria = [];
    end
    
    % Find dim and embedded dim
    ind_dim = [];
    ind_embedded = [];
    ndims = length(s.meta);
    for i = 1:ndims
        tmp = strcmp(s.meta{i}.dims_embedded,dim_name);
        if any(tmp)
            ind_dim = i;
            ind_embedded = find(tmp);
        end
    end
    
    if isempty(ind_dim) || isempty(ind_embedded)
        error('Chosen dimension not found');
    end
    
    vals = s.meta{ind_dim}.values_mat(:,ind_embedded);
    
    % Find chosen indices within dim
    if isempty(selection_criteria) || strcmp(selection_criteria,':')
        ind_selected = true(size(vals));
    else

        % Convert selection to cell if it's not already one, for simplicity
        if ~iscell(selection_criteria)
            selection_criteria = {selection_criteria};
        end

        % Identify selected
        ind_selected = false(size(vals));
        for i = 1:length(selection_criteria)
            sc = selection_criteria{i};
            if ~ischar(sc)          % If its not a string, assume we're just indexing
                ind_selected(sc) = true;
            else                    % If is a string, assume its a regular expression
                ind = [];
                func = @(x) isempty(string_matching_func(x,sc));
                ind = ~cellfun(func,vals);
                ind_selected(ind) = true;
            end
        end
    end
    
end