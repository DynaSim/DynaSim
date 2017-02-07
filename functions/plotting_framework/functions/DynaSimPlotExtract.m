

function [data_linear,ax,ax_names,time] = DynaSimPlotExtract (data)
    % Converts DynaSim structure to 1D cell array format. Later can use to
    % import to xPlt

    
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
    
    % ## Build a large linear list ##
    z=0;
    num_alllabels = 1;
    for i = 1:num_sims
        for k = 1:num_alllabels
            
            z=z+1;
            data_linear{z} = data(i).plot_files{1};
            
            % Number of parameter sweeps, plus populations, plus variables (Vm, state variables, functions, etc.)
            for j = 1:num_varied
                if isnumeric(data(1).(varied{j})); ax{j}(z) = params{j}(i);
                else
                    ax{j}{z} = params{j}{i};
                end
            end
        end
    end
    
    ax_names = varied;
    
end


