function xp = importLinearData(xp,X,varargin)
    % xp = importLinearData(X,axislabels1,...,axislabelsN)
    % xp = importLinearData(X,axislabels1,...,axislabelsN,'outputformat',Value)
    % Imports a linear array of data, converts it into a matrix
    % based on the supplied axislabels, and stores it in xp.data.
    % Also populates the xp.axis.values appropriately.
    % X - linear matrix or cell array containing input data
    % axislabels1 - linear matrix or cell array containing data labels for dimension1
    % ...
    % axislabelsN - linear matrix or cell array containing data labels for dimensionN
    % Value - Specifies storage format of the imported data, either 'cell' or 
    %         'numeric'. Note: if input data is a cell,
    %                output format must also be set to 'cell.'

    % Initialize
    axeslinear = varargin;

    % Set outputformat as needed based on data type of input format
    if iscell(X)
        outputformat = 'cell';          
    else
        outputformat = 'numeric';
    end

    % Check if final argument in varargin is a name/value pair
    if ischar(varargin{end-1})
        if strcmp(varargin{end-1},'format')
            outputformat = varargin{end};
            axeslinear = axeslinear(1:end-2);       % Remove the name-value pair from axeslinear
        end
    end


    % Error checking
    if ~isvector(X); error('X must be linear'); end
    N = length(X);
    Ndims = length(axeslinear);

    % Set up xp.axis
    for j = 1:Ndims
        xp.axis(j).values = unique(axeslinear{j});
        sz(j) = length(xp.axis(j).values);

        if isnumeric(axeslinear{j}(1))
            if any(isnan(axeslinear{j})) || any(isinf(axeslinear{j}))
                error('Axis cannot contain NaNs or Infs');
            end
        end
    end

    % Set up target matrix
    if strcmp(outputformat,'cell')
        xp.data={};
    else
        xp.data = zeros(sz);
    end

    % Set up xp.data -> Convert linear data into a multi dimensional matrix
    for i = 1:N
        % Get subscripts
        subs = cell(1,Ndims);
        for j = 1:Ndims
            if iscell(axeslinear{j})
                subs{j} = find(strcmp(axeslinear{j}{i},xp.axis(j).values));
            else
                subs{j} = find(axeslinear{j}(i) == xp.axis(j).values);
            end
        end

        % Add data to sparse cell array based on subscripts
        ind = sub2ind(sz,subs{:});
        xp.data(ind) = X(i);
    end


    % Increase size of cell array if necessary, to make it "square" so
    % that the subsequent reshape will work
    if strcmp(outputformat,'cell') && length(xp.data) < prod(sz)
        xp.data{prod(sz)} = [];
    end

    % Lastly convert the linear cell array to a matrix cell array
    xp.data = reshape(xp.data,sz);
end