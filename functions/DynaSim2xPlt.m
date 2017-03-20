
function xp = DynaSim2xPlt(data)
    % Convert DynaSim data structure to xp format

    data = CheckData(data);
    
    % Extract the data in a linear table format
    [data_table,column_titles,time] = Data2Table (data);

    % % Preview the contents of this table
    % %     Note: We cannot make this one big cell array since we want to allow
    % %     axis labels to be either strings or numerics.
    % previewTable(data_table,column_titles);

    % Import the linear data into an xPlt object
    xp = xPlt;
    X = data_table{1};                          % X holds the data that will populate the multidimensional array. Must be numeric or cell array.
    axislabels = data_table(2:end);             % Each entry in X has an associated set of axis labels, which will define its location in multidimensional space. **Must be numeric or cell array of chars only**
    xp = xp.importLinearData(X,axislabels{:});
    xp = xp.importAxisNames(column_titles(2:end));  % There should be 1 axis name for every axis, of type char.

    % Store metadata info
    meta = struct;
    meta.datainfo(1:2) = nDDictAxis;
    meta.datainfo(1).name = 'time(ms)';
    meta.datainfo(1).values = time;
    meta.datainfo(2).name = 'cells';
    meta.datainfo(2).values = [];
    meta.dynasim.labels = data.labels;
    meta.dynasim.model = data.model;
    meta.dynasim.simulator_options = data.simulator_options;
    meta.dynasim.time = data.time;
    meta.dynasim.varied = data.varied;
    xp.meta = meta;
    clear meta
end