function vars_pull(s,chosen_fields)
    % Pulls all fields in structure s into variables in local work space.
    % If chosen_fields is specified, only pulls in the chosen fields.
    if nargin < 2
        chosen_fields = fieldnames(s);
    end
    
    chosen_fields = chosen_fields(:)';

    for n = chosen_fields
        name = n{1};
        value = s.(name);
        assignin('caller',name,value);
    end
    
end