function xp = img2xPlt(data_img)
    % Load into DynaSim structure
    [data_table,column_titles] = ds.dataField2Table (data_img,'plot_files');

    % The entries in the first column contain the paths to the figure files.
    % There can be multiple figures associated with each simulation, which is
    % why these are cell arrays of strings.
    %disp(data_table{1}{1})
    %disp(data_table{1}{2})

    % Import the linear data into an xPlt object
    xp = xPlt;
    X = data_table{1}; axislabels = data_table(2:end);
    xp = xp.importLinearData(X, axislabels{:});
    xp = xp.importAxisNames(column_titles(2:end));

    % Squeeze out any empty dims that might have been introduced by the
    % above operations. This is necessary if xp was originally 2x1 (e.g.
    % varied1 x Dim 1) and then added populations and variables onto this
    % after Dim 1
    xp = xp.squeezeRegexp('Dim');
    
    % Set up metadata
    % Store metadata info (putting random info in here for time being)
    meta = struct;
    meta.datainfo(1:2) = nDDictAxis;
    meta.datainfo(1).name = 'time(ms)';
    meta.datainfo(1).values = 1:10;
    meta.datainfo(2).name = 'cells';
        cell_names = [1:5];
        cell_names_str = cellfunu(@(s) ['Cell ' num2str(s)], num2cell(cell_names));    
    if isfield(data_img(1),'varied')
        meta.dynasim.varied = data_img(1).varied;
    else
        % For case when nothing varied, insert some dummy data
        meta.dynasim.varied = {'Varied1'};
    end
    xp.meta = meta;
    
    
    
end
