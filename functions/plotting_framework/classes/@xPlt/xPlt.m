
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
        

        function [xp2, ro] = subset(xp,varargin)
            % Define variables and check that all dimensions are consistent
            % ro - if regular expressions are used, returns the index
            % values discovered by the regular expression.
            checkDims(xp);
            selection = varargin(:);
            Nd = ndims(xp);
            Na = length(selection);
            if Na ~= Nd
                error('Number of inputs must match dimensionality of xPlt.data');
            end
            
            % Convert selection to index if using regular expressions
            ro = {};
            for i = 1:Na
                if ischar(selection{i})
                    ro{i}(1,:) = xp.axis(i).values;
                    [selection{i} ro{i}(2,:)] = regex_lookup(xp.axis(i).values, selection{i});
                    
                end
            end
            
            % Initialize
            sz = size(xp);
            xp2 = xPlt;
            xp2.meta = xp.meta;
            
            % First update each axis with the corresponding selection and
            % convert empty cells to code for full range.
            for i = 1:Na

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
        
        function xp = importMeta(xp,meta_struct)
            xp.meta = meta_struct;
        end
        
        function [selection_out, startIndex] = findaxis(xp,str)
            % Returns the index of the axis with name matching str
            allnames = {xp.axis.name};
            [selection_out, startIndex] = regex_lookup(allnames, str);
        end
        
        function xp = mergeDims(xp,dims2pack)
            
            xp.checkDims;
            Nd2p = length(dims2pack);
            sz = size(xp.data);
            Nmerged = prod(sz(dims2pack));       % Final number entries in the merged dimension
            
            
            % % First, do axis names
            
            % Get cell array of all linearized axis values.
            inds = 1:Nmerged;
            [subs{1:Nd2p}] = ind2sub(sz(dims2pack),inds);
            temp = cell(1,Nd2p);
            for i = 1:Nd2p
                for j = 1:Nmerged
                    ax = xp.axis(dims2pack(i));
                    currval = ax.values(subs{i}(j));
                    temp{i}(j) = currval;
                end
            end
            
            % Compress these into a single string to be used as new value name
            tempstr = {};
            for i = 1:Nd2p
                for j = 1:Nmerged
                    if iscell(temp{i}); currval = temp{i}{j};
                    else
                        currval = num2str(temp{i}(j));  % If it's not a string, convert it to one.
                    end
                    if i == 1; tempstr{j} = currval;
                    else tempstr{j} = [tempstr{j} '_' currval];
                    end
                end
            end
            
            % Finally drop this into the values entry for the new "merged"
            % axis
            xp.axis(dims2pack(1)).values = tempstr;
            xp.axis(dims2pack(1)).astruct.premerged_values = temp;
            
            
            % Give it a new axis name, reflecting the merger of all the
            % others
            allnames = {xp.axis(dims2pack).name};
            allnames = cat(1,allnames(:)',repmat({'_'},1,length(allnames)));
            xp.axis(dims2pack(1)).name = strcat(allnames{1:end-1});
            
            % Clear the remaining axes names
            for i = 2:Nd2p
                xp.axis(dims2pack(i)).name = ['Dim ' num2str(dims2pack(i))];
                xp.axis(dims2pack(i)).values = 1;
            end


            % % Now, work on xp.data
            dat = xp.data;
            xp.data = [];       % Clear xp.data for now, to save memory.
            
            % Figure out which dimensions were *NOT* Targeted for the merge
            Nd = ndims(dat);
            alldims = 1:Nd;
            ind_chosen = false(size(alldims));
            for i = 1:length(dims2pack)
                ind_chosen = ind_chosen | alldims == dims2pack(i);
            end
            ind_unchosen = ~ind_chosen;
            dims_remaining = find(ind_unchosen);
            
            % Bring the dims to be merged to the front
            dat = permute(dat,[dims2pack,dims_remaining]);
            
            % REshape these into a single dim
            sz = size(dat);
            dat = reshape(dat,[ prod(sz(1:Nd2p)), ones(1,Nd2p-1), sz(Nd2p+1:end) ]);
            
            % Undo the earlier permute, and put back into xp.data
            dat = ipermute(dat,[dims2pack,dims_remaining]);
            xp.data = dat;
        end
        
        function xp = packDim(xp,dim_src,dim_target)
            
            if nargin < 3; dim_target = dim_src; end
            checkDims(xp);
            
            % Make sure that xp.data is a cell array
            if ~iscell(xp.data); error('xPlt.data must be a cell array.'); end
            
            % Make sure that xp.data is a matrix
            temp = cellfun(@ismatrix,xp.data);
            if any(temp(:) ~= 1); error('xPlt.data must contain only matrices'); end
            % % To do: implement this so it works with cell arrays and xPlt
            % classes in the future too
            
            % Make sure target dimension in xPlt.data is a singleton
            temp = cellfun(@(x) size(x,dim_target),xp.data);
            if any(temp(:) ~= 1); error('Target dimension in xPlt.data needs to be size 1. Try reshaping contents of xPlt.data or choosing a different target dimension.'); end
            clear sz_targets
            
            % Bring chosen dimension to the front. Thus, we will be
            % merging along rows.
            Nd = ndims(xp.data);
            xp.data = permute(xp.data,[dim_src, 1:dim_src-1, dim_src+1:Nd]);
            
            % Temporarily linearize all other dimensions.
            sz0 = size(xp.data);
            xp.data = reshape(xp.data,sz0(1),prod(sz0(2:end)));
            
            % Add NaNs where needed
                % Note: to understand what this is doing, it really helps
                % to draw a diagram!
            sz = size(xp.data);
            empties = cellfun(@isempty,xp.data);    % 2D matrix with 1's marking empty cells
            se = sum(empties,1);                    % Number of empties per column in this matrix
            bad_inds = se ~= 0 & se ~= sz(2);     % Good columns must be either all empty or all non-empty
            
            if any(bad_inds)
                fprintf('Note: Empty entries found along collapsing dim. Using NaNs as placeholders to fill out the matrix. \n');
                bi = find(bad_inds);
                for j = 1:length(bi)                    % Sweep through bad columns
                    curr_bad = find(empties(:,j));      % Empties along this particular column
                    curr_good = find(~empties(:,j));    % Non-empties along this particular column.
                    for i = 1:length(curr_bad)
                        % Populate the empty cells with matrices of NaNs
                        % that are the same dimensionality as the first
                        % good entry.
                        xp.data{curr_bad(i),bi(j)} = NaN*ones(size(xp.data{curr_good(1),bi(j)}));
                    end
                end
            end
            
            % Check that sizes and dimensionalities are compatible
            data_ndims = cellfun(@ndims,xp.data,'UniformOutput',true);
            if any(any(data_ndims ~= repmat(data_ndims(1,:),sz(1),1),1),2)
                error('Dimensions of xPlt.data not uniform along packing dimensions.');
            end
           
            
            data_sz = cellfun(@size,xp.data,'UniformOutput',false);
            data_sz_firsts = repmat(data_sz(1,:),sz(1),1);
            myfunc = @(x,y) any(x(:) ~= y(:));
            bool_size_mismatch = cellfun(myfunc,data_sz,data_sz_firsts);
            if any(bool_size_mismatch(:))
                error('Sizes of xPlt.data are not uniform along packing dimension. (This usually results form trying to combine populations with different numbers of cells.');
            end
            
            % Now, pack the data.
%             for j = 1:sz(2)
%                 % permute so that contents of each cell has the selected
%                 % dimenison along dimension 1
%                 for i = 1:sz(1)
%                     dat_curr = xp.data{i,j};
%                     Ndimdat = ndims(dat_curr);
%                     dat_curr = permute(dat_curr,[dim, 1:dim-1, dim+1:Ndimdat]); %whos dat_curr
%                     % dat_curr = permute(dat_curr,[2:dim, 1, dim+1:Ndimdat]);    %whos dat_curr       % Uncomment this stuff to make sure it permutes back
%                     xp.data{i,j} = dat_curr;
%                 end
%                 
%                 % Concatenate the data along dimension 1
%                 data_new{j} = cat(1,xp.data{j,:});
%                 
%                 % Permute the contents of each cell back to original
%                 for i = 1:sz(1)
%                     dat_curr = data_new{j};
%                     dat_curr = permute(dat_curr,[2:dim, 1, dim+1:Ndimdat]);     %whos dat_curr
%                     data_new{j} = dat_curr;
%                 end
%                 
%             end
        
            for j = 1:sz(2)
                xp.data{1,j} = cat(dim_target,xp.data{:,j});
                xp.data(2:end,j) = cell(sz(1)-1,1);     % Set remainder to empty
            end
            
            xp.data = xp.data(1,:);         % Keep only 1st dimension;
            sz(1) = 1;
            sz0(1) = 1;

            % Lastly, restore original dimensions
            % of xPlt.data
            xp.data = reshape(xp.data,sz0);
            xp.data = permute(xp.data,[2:dim_src, 1, dim_src+1:Nd]);
            
            % Also, update xp.axis
%             xp.axis(dim_src).values = 1;
            xp = setAxisDefaultNames(xp,dim_src);
            xp.axis(dim_src).name = ['Dim ' num2str(dim_src)];
            
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

            % Trim away excess values in axes
            if Na > Nd
                for i = Nd+1:Na
                    if length(xp.axis(i).values) > 1
                        xp.axis(i).values = xp.axis(i).values(1);
                        fprintf(['Extra values found in axis #' num2str(i) ' ' xp.axis(i).name '. Trimming \n']);
                    end
                end
                
            end
            
        end
        
        
        function checkDims(xp)

            % Note, fixAxes fixes everything automatically.
            % Only call checkDims if you want to
            % be alerted to mismatches, but not to correct them. Use fixAxes to
            % automatically correct everything.
            
            sz = size(xp);
            Nd = length(sz);
            Na = length(xp.axis);



            if Nd ~= Na
                error('checkDims: Error found! Number of dimensions in xPlt.data does not equal number of axes');
            end

            for i = 1:Na
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

        
        % % % % % % % % % % % OVERLOADED FUNCTIONS % % % % % % % % % % %
        
        function varargout = size(xp,varargin)
            % Overrides normal size command.
%             for j = 1:length(xp.axis)
%                 sz(j) = length(xp.axis(j).values);
%             end

            [varargout{1:nargout}] = size(xp.data,varargin{:});
            
            % Add singleton dimensions as needed
            if nargout <= 1 && nargin == 1
                sz = varargout{1};
                Nd = ndims(xp.data);
                Na = length(xp.axis);
                if Nd < Na
                    sz_axis = cellfun(@length,{xp.axis.values});
                    if any(sz_axis(Nd+1:Na) > 1); error('Dimension mismatch between size(xp.data) and length(xp.axis.values). Run getaxisinfo or fixAxes.'); end
                    sz(Nd+1:Na) = 1;
                end
                varargout{1} = sz;
            end
            
        end
        
        function Nd = ndims(xp)
            Nd = ndims(xp.data);
            
            % If there are more axes than Nd, this could be because there
            % are a bunch of singleton axes. If there are, then we can
            % update this returned value.
            Na = length(xp.axis);
            if Nd < Na
                sz_axis = cellfun(@length,{xp.axis.values});
                if any(sz_axis(Nd+1:Na) > 1); error('Dimension mismatch between size(xp.data) and length(xp.axis.values). Run getaxisinfo or fixAxes.'); end
                Nd = Na;
            end

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
            
            % If data is bigger than a matrix, squeeze out dimensions that
            % are of size 1.
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
            else
                % Otherwise, if data is a matrix, remove all axis beyond
                % the first two. These should only be size 1 (e.g. "name"
                % axes anyways)
%                 szA = cellfun(@length,{xp.axis.values});
%                 ind = szA~=1;
%                 ind(1:2) = true;
                xp.axis = xp.axis(1:2);
            end
            
            % Make sure everything is good before returning.
            xp = xp.fixAxes;
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




function xp = setAxisDefaultNames(xp,dim)
    % Sets xp.axis(i) to default values
    
    % Get desired size of dataset
    sz_dim = size(xp.data,dim);
    
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
        ax_curr.values = 1:sz_dim;                                                   % Populate with numerics
    else
        % Otherwise, make sure dimensionality is correct. If not, update it
        % missing entries with default names.
        N = length(ax_curr.values);
        if N < sz_dim
            if isnumeric(ax_curr.values)
                for j = N:sz_dim; ax_curr.values(j) = j; end
            else
                for j = N:sz_dim; ax_curr.values{j} = num2str(j); end
            end
        end
        
        if N > sz_dim
            %ax_curr.values = ax_curr.values(1:sz(dim));
            ax_curr.values = 1:sz_dim;                                                   % Populate with genetic numerics
        end
    end
    
    xp.axis(dim) = ax_curr;
end


function [selection_out, startIndex] = regex_lookup(vals, selection)
    if ~ischar(vals{1}); error('xPlt.axis.values must be strings when using regular expressions');
    end
    if ~ischar(selection); error('Selection must be string when using regexp');
    end
    
    startIndex = regexp(vals,selection);
    selection_out = logical(~cellfun(@isempty,startIndex));
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
