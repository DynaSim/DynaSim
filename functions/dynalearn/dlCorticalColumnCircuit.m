function y = dlCorticalColumnCircuit(name, layers, celltypes, counts, connections, noises)

    y = struct();

    if ~exist('name', 'var')

        disp("Warning! Name not specified.");
        y.name = "Network0";

    else

        y.name = name;

    end

    if ~exist('layers', 'var')

        disp("Warning! Number of layers not specified.");
        y.layers = 3;

    else

        y.layers = layers;

    end

    if ~exist('celltypes', 'var')

        disp("Cell types not specified; Default : [""E"", ""PV"", ""CB"", ""CR""]");
        y.celltypes = ["E", "PV", "CB", "CR"];
        y.celltypesCount = length(y.celltypes);

    else

        y.celltypes = string(celltypes);
        y.celltypesCount = length(y.celltypes);

    end

    if ~exist('counts', 'var')

        disp("Warning! Size of the populations (Layer x Celltype) not specified.")
        y.counts = ones(y.layers, length(y.celltypes));

    else

        y.counts = counts;

    end

    if ~exist('connections', 'var')

        disp("Warning! connections between populations not specified.");
        y.connections = zeros(y.layers*length(y.celltypes), y.layers*length(y.celltypes));

    else

        y.connections = connections;

    end

    if ~exist('noises', 'var')

        y.noises = ones(1, y.celltypesCount)*20;

    else

        y.noises = noises;

    end

end