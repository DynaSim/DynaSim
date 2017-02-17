

function [xp] = DynaSim2xPlt (data)
    % Converts DynaSim structure to cell array format for use in xPlt
    % framework. This function directly copies the data into the cell,
    % rather than linearizing and using importLinearData.
    
    
    xp = xPlt;
    curr_dim = 0;
    

    % Add Time variable
    xp.meta.time = data(1).time;
    
    
    
    % ## VARIED parameter sweeps ##
    varied=data(1).varied;
    num_varied=length(varied); % number of model components varied across simulations
    num_sims=length(data); % number of data sets (one per simulation)
    param_mat=zeros(num_sims,num_varied); % values for each simulation
    param_cell=cell(1,num_varied); % unique values for each parameter
    
    for j=1:num_varied
        
        if isnumeric(data(1).(varied{j}))
            param_mat(:,j)=[data.(varied{j})]; % values for each simulation
        else
            for i = 1:length(data)
                param_mat{i,j} = data(i).(varied{j});    %vals for each sim
            end
        end
        param_cell{j}=unique([data.(varied{j})]); % unique values for each parameter
        
        fld = varied{j};
        for k = 1:length(param_cell{j})
            str{k} = [fld '=' num2str(param_cell{j}(k))];
        end
        
        % ## Store parameter sweep in a dimension ##
        curr_dim = curr_dim + 1;
        xp.axis(curr_dim) = xPltAxis;
        xp.axis(curr_dim).name = varied{j};
        xp.axis(curr_dim).values = param_cell{j};
    end
    

    % ## FIELDS Get all fields of data ##
    % Get metadata for all data fields (e.g. populations / currents / state
    % variables) 
    labels = data(1).labels;
    labels = labels(cellfun(@isempty,strfind(labels,'time')));      % Remove time from labels
    num_alllabels = length(labels);
    
    % Determine all unique populations
    pop_names={data(1).model.specification.populations.name}; % list of populations
    
    % Build list populations and variables
    ind = strfind(labels,'_');
    func1 = @(x,y) x(1:y(1)-1);
    func2 = @(x,y) x(y(1)+1:end);
    pops = cellfun(func1,labels,ind,'UniformOutput',0);
    vars = cellfun(func2,labels,ind,'UniformOutput',0);
    
    % Determine all unique variables
    var_names = unique(vars);
    
    % ## Store populations in a dimension ##
    curr_dim = curr_dim + 1;
    xp.axis(curr_dim) = xPltAxis;
    xp.axis(curr_dim).name = 'populations';
    xp.axis(curr_dim).values = pop_names;

    
    % ## Store variables in a dimension ##
    curr_dim = curr_dim + 1;
    xp.axis(curr_dim) = xPltAxis;
    xp.axis(curr_dim).name = 'variables';
    xp.axis(curr_dim).values = var_names;

    
    % Calculate sizes of each dimension
    sz = xp.get_sz;
    
    % Pack in to data structure
    inds = zeros(1,length(xp.axis));
    z=0;
    for i = 1:num_sims
        for k = 1:num_alllabels
            z=z+1;
            
            % Get dimensions for varied parameter sweeps
            curr_dim = 0;
            for j = 1:num_varied
                curr_dim = curr_dim + 1;
                if isnumeric(param_mat(1))
                    inds(curr_dim) = find(xp.axis(curr_dim).values == param_mat(i,j));
                else
                    inds(curr_dim) = strcmp(xp.axis(curr_dim).values,param_mat(i,j));
                end
            end
            
            % Get dimensions for populations
            curr_dim = curr_dim + 1;
            inds(curr_dim) = find(strcmp(xp.axis(curr_dim).values,pops{k}));
            
            % Get dimensions for variables
            curr_dim = curr_dim + 1;
            inds(curr_dim) = find(strcmp(xp.axis(curr_dim).values,vars{k}));
            
            % Get the linear index
            indsc = num2cell(inds);
            linear_ind(z) = sub2ind(sz,indsc{:});
            
            % Store the data
            xp.data{linear_ind(z)} = data(i).(labels{k});
            
%             warning('Re-write this using customized import commands');
%             warning('Write a check function to compare size of xp.axis data to size of xp.data');
%             warning('Also write functions to generate other entries besides sz, such as list of populations, dimensions, etc - or to return these in pleasing formats');
%             warning('Write a squeeze function');
%             warning('Split this function up into several sub functions that return all the meta data gathered here!');
        end
    end
    
    % Increase size of data matrix if necessary, to make it "square" so
    % that the subsequent reshape will work
    if length(xp.data) < prod(sz)
        xp.data{prod(sz)} = [];
    end

    % Finally reshape the data matrix to the appropriate dimensions
    xp.data = reshape(xp.data,sz);
    
end


