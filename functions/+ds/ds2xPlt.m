function xp = ds2xPlt(data, varargin)
    % Convert DynaSim data structure to xp format

    data = ds.checkData(data, varargin{:});
    
    % Extract the data in a linear table format
    [data_table,column_titles,time] = ds.data2Table (data);

    % % Preview the contents of this table
    % %     Note: We cannot make this one big cell array since we want to allow
    % %     axis labels to be either strings or numerics.
    % ds.previewTable(data_table,column_titles);

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
        cell_names = [1:max(cellfun(@(x) size(x,2),xp.data(:)))];
        cell_names_str = cellfunu(@(s) ['Cell ' num2str(s)], num2cell(cell_names));
    meta.datainfo(2).values = cell_names_str;
    meta.dynasim.labels = data(1).labels;
    meta.dynasim.model = data(1).model;
    meta.dynasim.simulator_options = data(1).simulator_options;
    meta.dynasim.time = data(1).time;
    if isfield(data(1),'varied')
        meta.dynasim.varied = data(1).varied;
    else
        % For case when nothing varied, insert some dummy data
        meta.dynasim.varied = {'Varied1'};
    end
    xp.meta = meta;
    clear meta
end
