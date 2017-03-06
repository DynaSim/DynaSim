function obj = importLinearData(obj,X,varargin)
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


% %     Don't need this option for now.
%     % Check if final argument in varargin is a name/value pair
%     if ischar(varargin{end-1})
%         if strcmp(varargin{end-1},'format')
%             outputformat = varargin{end};
%             axeslinear = axeslinear(1:end-2);       % Remove the name-value pair from axeslinear
%         end
%     end



    % Error checking - X must be linear
    if ~isvector(X); error('X must be linear'); end


    % Error checking - X must be cell or numeric
    obj2 = obj.reset;
    outputformat = obj2.validateInputs(X, 'data');


    N = length(X);

    % Error checking - each entry in axislinear must be either numeric or
    % cell. If it's a cell, all entries must char.
    temp = obj2.validateInputs(axeslinear,'axis');


    Ndims = length(axeslinear);

    % Set up xp.axis
    for j = 1:Ndims
        obj.axis(j).values = unique(axeslinear{j});
        sz(j) = length(obj.axis(j).values);

        if isnumeric(axeslinear{j}(1))
            if any(isnan(axeslinear{j})) || any(isinf(axeslinear{j}))
                error('Axis cannot contain NaNs or Infs');
            end
        end
    end

    % Set up target matrix
    switch outputformat
        case 'cell'
            obj.data=cell(sz);
%         case 'string'
%             xp.data = repmat(string(''),sz);
        case 'numeric'
            obj.data = zeros(sz);
        otherwise
            error('Case not implemented');
    end

    % Set up xp.data -> Convert linear data into a multi dimensional matrix
    for i = 1:N
        % Get subscripts
        subs = cell(1,Ndims);
        for j = 1:Ndims
            if iscellstr(axeslinear{j})
                subs{j} = find(strcmp(axeslinear{j}{i},obj.axis(j).values));
            else
                subs{j} = find(axeslinear{j}(i) == obj.axis(j).values);
            end
        end

        % Add data to sparse cell array or matrix based on subscripts
            % Note: Need to find a good way for dealing with duplicate
            % rows. Right now, default behavior is to overwrite
        obj.data(subs{:}) = X(i);

    end

end


