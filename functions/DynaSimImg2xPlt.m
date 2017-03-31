function xp = DynaSimImg2xPlt(data_img)
    % Load into DynaSim structure
    [data_table,column_titles] = DataField2Table (data_img,'plot_files');

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
    
    % Add dummy population and variable dimensions since the code below
    % expects it
    xd = xp.data;
    xv = xp.exportAxisVals;
    xn = xp.exportAxisNames;
    
    Nd = ndims(xp.data);
    if Nd == 2 && size(xp.data,2) == 1      % If it's a 2x1 matrix, really this is a vector not a matrix, so set Nd to 1.
        Nd = 1;
    end
    
    xv(Nd+1:Nd+2) = {{'Pop1'},{'X'}};
    xn(Nd+1:Nd+2) = {'populations','variables'};
    
    xp = xp.importData(xd,xv);
    xp = xp.importAxisNames(xn);
    clear xd xv xn
    
    % Set up metadata
    % Store metadata info
    meta = struct;
    meta.datainfo(1:2) = nDDictAxis;
    meta.datainfo(1).name = 'time(ms)';
    meta.datainfo(1).values = 1:10;
    meta.datainfo(2).name = 'cells';
        cell_names = [1:5];
        cell_names_str = cellfunu(@(s) ['Cell ' num2str(s)], num2cell(cell_names));
    xp.meta = meta;
end