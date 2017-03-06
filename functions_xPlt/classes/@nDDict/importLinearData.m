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
    N = length(X);
    
    % Error checking - X must be cell or numeric
    if iscell(X)
        outputformat = 'cell';          
    elseif isnumeric(X)
        outputformat = 'numeric';
%     elseif isstring(X)
%         outputformat = 'string';
    else
        error('Input X must be cell or numeric');
    end

    % Error checking - each entry in axislinear must be either numeric or
    % cell. If it's a cell, all entries must char.
    for i = 1:length(axeslinear)
        out = validate_axis(axeslinear{i});
        if strcmp(out,'unknown')
            error(['axislabels' num2str(i) ' must be of type numeric or cell array of character vectors']);
        end
    end
    
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
    switch outputformat
        case 'cell'
            xp.data=cell(sz);
%         case 'string'
%             xp.data = repmat(string(''),sz);
        case 'numeric'
            xp.data = zeros(sz);
        otherwise
            error('Case not implemented');
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
        
        % Add data to sparse cell array or matrix based on subscripts
            % Note: Need to find a good way for dealing with duplicate
            % rows. Right now, default behavior is to overwrite
        xp.data(subs{:}) = X(i);
        
    end

end



function ax_type = validate_axis(axlinear)

    if isnumeric(axlinear)
        ax_type = 'numeric';
%     elseif isstring(axlinear)
%         ax_type = 'string';
    elseif iscell(axlinear)
        if is_cellarray_of_chars(axlinear)
            ax_type = 'cell_of_chars';
        else
            ax_type = 'unknown';
        end
    else
        ax_type = 'unknown';
    end
    

    

end

function out = is_cellarray_of_chars(c)
    out = all(cellfun(@ischar,c));
end

% This code is for enabling the use of string  data types as well.
% Disabling for now, since string only works in 2016b and backwards
% compatibility is complicated.
% function varargout = string(varargin)
%     % Overload string function for backwards compatibility
%     if exist('string')
%         [varargout{1:nargout}] = builtin('string',varargin{:});
%     else
%         c = varargin{1};
%         if iscell(c)
%             is_cellarray_of_chars(c);
%             out = c;
%         elseif ischar(c)
%             out = {c};
%         else
%             error('Unknown input type');
%         end
%         
%         varargout{1} = out;
%     end
% end
% 
% function varargout = isstring(varargin)
%     if exist('isstring','builtin')
%         [varargout{1:nargout}] = builtin('isstring',varargin{:});
%     else
%         c = varargin{1};
%         if iscell(c)
%             out = is_cellarray_of_chars(c);
%         elseif ischar(c)
%             out = true;
%         else
%             out=false;
%         end
%         varargout{1} = out;
%     end
% end
