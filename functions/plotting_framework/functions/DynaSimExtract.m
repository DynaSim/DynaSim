

function [data_linear,ax,ax_names,time] = DynaSimExtract (data)
    % Converts DynaSim structure to 1D cell array format. Later can use to
    % import to xPlt

    
    % Extract Time variable
    time = data(1).time;
    
    % ## VARIED parameter sweeps ##
    varied=data(1).varied;
    num_varied=length(varied); % number of model components varied across simulations
    num_sims=length(data); % number of data sets (one per simulation)
    %param_mat=zeros(num_sims,num_varied); % values for each simulation
    
    for j=1:num_varied
        if isnumeric(data(1).(varied{j}))
            params{j} = [data.(varied{j})]; % values for each simulation
        else
            for i = 1:length(data)
                params{j}{i} = data(i).(varied{j});    %vals for each sim
            end
        end
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
    
    % ## Build a large linear list ##
    z=0;
    for i = 1:num_sims
        for k = 1:num_alllabels
            
            z=z+1;
            data_linear{z} = data(i).(labels{k});
            
            % Number of parameter sweeps, plus populations, plus variables (Vm, state variables, functions, etc.)
            for j = 1:num_varied
                if isnumeric(data(1).(varied{j})); ax{j}(z) = params{j}(i);
                else
                    ax{j}{z} = params{j}{i};
                end
            end
            
            ax{num_varied+1}{z} = pops{k};
            
            j=j+1;
            ax{num_varied+2}{z} = vars{k};
            
        end
    end
    
    ax_names = varied;
    ax_names{num_varied+1} = 'populations';
    ax_names{num_varied+2} = 'variables';
end


