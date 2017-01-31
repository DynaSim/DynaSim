
classdef xPlt

    properties
        data               % Storing the actual data (multi-dimensional matrix or cell array)
        axis = xPltAxis    % 1xNdims - array of xPltAxis classes for each axis. Ndims = ndims(data)
        meta = struct;     % Metadata about stuff that's stored in data
    end
    
    
    methods

        function xp = xPlt(data,axis,meta)
            if exist('data','var')
                xp.data = data;
            end
            
            if exist('axis','var')
                xp.axis = axis;
            end
            
            if exist('meta','var')
                xp.meta = meta;
            end
            
        end
        
        
        
        function xp2 = subset(xp,varargin)
            % To do: Need to make this work with ways to select the subset
            % based on value name, rather than just by index.
            
            % Define variables and check that all dimensions are consistent
            checkDims(xp);
            selection = varargin(:);
            Nd = ndims(xp.data);
            N = length(selection);
            if N ~= Nd
                error('Number of inputs must match dimensionality of xPlt.data');
            end
            
            % Convert selection to index if using regular expressions
            for i = 1:length(selection)
                if ischar(selection{i})
                    selection{i} = regex_lookup(xp.axis(i).values, selection{i});
                end
            end
            
            % Initialize
            sz = size(xp);
            xp2 = xPlt;
            xp2.meta = xp.meta;
            
            % First update each axis with the corresponding selection and
            % convert empty cells to code for full range.
            for i = 1:length(selection)

                if isempty(selection{i})
                    selection{i} = 1:sz(i);
                end
                
                xp2.axis(i) = xp.axis(i);       % Import axis information
                xp2.axis(i).values = xp.axis(i).values(selection{i});   % Overwrite values field; leave everything else the same.
            end
            
            % Update the data
            xp2.data = xp.data(selection{:});
            
            % Corrects number of axes. The above code automatically
            % converts xp.data from MxNx1x1 to MxN, whereas axis will stay
            % as MxNx1x1 (e.g. length of 4). Thus, fixAxes corrects this.
            xp2 = fixAxes(xp2);
            
        end
        
        function xp = importAxisNames(xp,ax_names)
            Nd = ndims(xp.data);
            if nargin < 2
                ax_names = cellfun(@num2str,num2cell(1:5),'UniformOutput',0);
            end
            
            if length(ax_names) ~= Nd
                error('Mismatch between number of axis names supplied and number of dimensions in dataset'); end

            for i = 1:ndims(xp.data)
                xp.axis(i).name = ax_names{i};
            end
        end
        
        function xp = packdims(xp,dims2pack)
            error('this is incomplete');
            % Calculate dims
            Nd = ndims(xp.data);
            alldims = 1:Nd;
            ind_chosen = false(size(alldims));
            for i = 1:length(dims2pack)
                ind_chosen = ind_chosen | alldims == dims2pack(i);
            end
            ind_unchosen = ~ind_chosen;
            dims_remaining = find(ind_unchosen);
            
            xp = xp.permute([dims2pack,dims_remaining]);
            
            sz = size(xp);
            xp.data = reshape(xp.data,[]);
        end
        
        
        
        function out = getaxisinfo(xp)
            % If no output arguments, prints axis info to the screen. If
            % output arguments are supplied, returns this information as a
            % strin.g
            
            if nargout > 0
                out = '';
            end
            
            for i = 1:length(xp.axis)
                
                out1 = xp.axis(i).getaxisinfo;
                temp = '';
                
                if nargout > 0
                    out = [out, temp, out1, '; ' ];
                else
                    temp = ['Axis ', num2str(i), ': '];
                    fprintf([temp, out1, '\n']);
                end
            end
            
            % Lastly output a summary of dimensionality comparing xPlt.axis
            % and xPlt.data. These should match up.
            if nargout == 0
                fprintf(['xPlt.axis dimensionality ' num2str(cellfun(@length,{xp.axis.values})) '\n']);
                fprintf(['xPlt.data dimensionality ' num2str(size(xp.data)) '\n']);
            end
        end
        
        % % % % % % % % % % % HOUSEKEEPING FUNCTIONS % % % % % % % % % % %
        function xp = fixAxes(xp)
            % This function forces the xPlt axis data to be updated to
            % match the dimensions of the data structure.
            % The convention of xPlt is to have always follow MATLAB
            % conventions for dimensionality. Thus, the size(xp.data)
            % command is used to determine the correct number of axes, and
            % axis is adjusted to match, adding or removing dimensions as
            % needed. If you are getting errors when running checkDims, you
            % should run this command.
            % 
            % Since xPlt.data and xPlt.axis can be manually edited by the
            % user, they can become mismatched in terms of their
            % dimensionality, this command can be used after making such
            % manual edits.
            %
            % Additionally, through ordinary useage of certain xPlt
            % functions, such as subset and squeeze, dimensions can get
            % mismatched. This happens because of certain conventions
            % MATLAB follows for treating matrices and cell arrays. For
            % example, MATLAB adds a trailing dimension onto column vectors
            % (e.g. the dimensionality is 2x1 instead of just 2).
            % Additionally, it concatenates off additional dimensions
            % beyond two, so a 2x1x1 is reduced to 2x1 automatically,
            % without need for a squeeze. This command forces xPlt.axis to
            % follow these conventions.

            Nd = ndims(xp.data);    
            Na = length(xp.axis);

            % Sweep through all axes and make sure dimensions are correct.
            % Add new axes if needed, up to Nd.
            for i = 1:Nd
                xp = setAxisDefaultNames(xp,i);  % Sets axis #i to the default name
            end

            % Trim away excess axes
            if Na > Nd
                xp.axis = xp.axis(1:Nd);
            end
            
        end

        
        % % % % % % % % % % % OVERLOADED FUNCTIONS % % % % % % % % % % %
        
        function varargout = size(xp,varargin)
            % Overrides normal size command.
%             for j = 1:length(xp.axis)
%                 sz(j) = length(xp.axis(j).values);
%             end

            [varargout{1:nargout}] = size(xp.data,varargin{:});
            
        end
        
        function xp = permute(xp,order)
            xp.data = permute(xp.data,order);
            xp.axis = xp.axis(order);
        end
        
        
        function xp = transpose(xp)
            checkDims(xp);
            Nd = ndims(xp.data);
            
            if Nd > 2; error('Can only transpose data with at most 2 dimensions');
            end
            
            xp.data = xp.data';
            xp.axis([1,2]) = xp.axis([2,1]);        % Axis should always be at least length=2.
        end
        
        function xp = squeeze(xp)
            % This is just like MATLAB's normal squeeze command. However,
            % there is one key difference:
            % Normally, if squeeze operates on a 1xN matrix, it will leave
            % it as 1xN. This function forces it to always return as Nx1
            
            checkDims(xp);
            
            % If it's bigger than a matrix, squeeze out dimensions that are
            % of size 1.
            sz = size(xp.data);
            if length(sz) > 2
                ind = sz~=1;
                xp.axis = xp.axis(ind);

                % Now squeeze xp.data
                xp.data = squeeze(xp.data);         % Normal squeeze command

%                 % Lastly, if the result is a row vector, force it to be a
%                 % column vector
%                 if isvector(xp.data) && ~iscolumn(xp.data)
%                     xp.data = xp.data';
%                 end
            end
            
            % Make sure everything is good before returning.
            checkDims(xp);
        end
        
        % % % % % % % % % % % END % % % % % % % % % % %
    end
end


function output = inheritObj(output,input)
    % Merges contents of input into output.
    C = metaclass(input);
    P = C.Properties;
    for k = 1:length(P)
        if ~P{k}.Dependent
            output.(P{k}.Name) = input.(P{k}.Name);
        end
    end
end




function xp = checkDims(xp)

    % Note, fixAxes fixes everything automatically now, thus rendering the
    % rest of this function pointless. Only call checkDims if you want to
    % be alerted to mismatches, but not to correct them. Use fixAxes to
    % automatically correct everything in the background.
    %xp = fixAxes(xp);

    sz = size(xp.data);
    Nd = length(sz);
    N = length(xp.axis);
    
    
    
    if Nd ~= N
        error('checkDims: Error found! Number of dimensions in xPlt.data does not equal number of axes');
    end

    for i = 1:N
        Nvalues_in_axis = length(xp.axis(i).values);
        if Nvalues_in_axis ~= sz(i)
            fprintf(['checkDims: Error found! Size of dimension ',num2str(i), ...
                ' is ', num2str(sz(i)) , ...
                '. But corresponding axis \"', xp.axis(i).name , ...
                '\" has ', num2str(Nvalues_in_axis) , ' elements. Dimension mismatch. \nTry running xPlt.getaxisinfo and then xPlt.fixAxes. \n' ]);
            error(' ');
        end
    end
end

function xp = setAxisDefaultNames(xp,dim)
    % Sets xp.axis(i) to default values
    
    % Get desired size of dataset
    sz = size(xp.data);
    
    % If axis doesn't already exist, create it. Otherwise, copy existing.
    if length(xp.axis) < dim
        ax_curr = xPltAxis;
    else
        ax_curr = xp.axis(dim);
    end
    
    % Name it if necessary
    if isempty(ax_curr.name)
        ax_curr.name = ['Dim ' num2str(dim)];
    end
    
    % If values is empty, add default values.
    if isempty(ax_curr.values)
        %ax_curr.values = cellfun(@num2str,num2cell(1:sz(i)),'UniformOutput',0);     % Populate with strings
        ax_curr.values = 1:sz(dim);                                                   % Populate with numerics
    else
        % Otherwise, make sure dimensionality is correct. If not, update it
        % missing entries with default names.
        N = length(ax_curr.values);
        if N < sz(dim)
            if isnumeric(ax_curr.values)
                for j = N:sz(dim); ax_curr.values(j) = j; end
            else
                for j = N:sz(dim); ax_curr.values{j} = num2str(j); end
            end
        end
        
        if N > sz(dim)
            %ax_curr.values = ax_curr.values(1:sz(dim));
            ax_curr.values = 1:sz(dim);                                                   % Populate with genetic numerics
        end
    end
    
    xp.axis(dim) = ax_curr;
end


function selection_out = regex_lookup(vals, selection)
    if ~ischar(vals{1}); error('xPlt.axis.values must be strings when using regular expressions');
    end
    if ~ischar(selection); error('Selection must be string when using regex');
    end
    
    selection_out = regexp(vals,selection);
    selection_out = logical(~cellfun(@isempty,selection_out));
    selection_out = find(selection_out);
    
end

% function varargout = size2(varargin)
%     [varargout{1:nargout}] = size(varargin{:});
%     if nargout == 1
%         sz = varargout{1};
%         if length(sz) == 2 && sz(2) == 1
%             sz = sz(1);
%         end
%         varargout{1} = sz;
%     end
% end
