function model = dlConnectModels(models, connections)

    m = length(connections);
    
    model = dsCombineSpecifications(models{:});

    for i = 1:m

        c = length(model.connections) + 1;
        model.connections(c) = connections{1, i};

    end

    fprintf("\n->Models Connected and Merged.\n");

end