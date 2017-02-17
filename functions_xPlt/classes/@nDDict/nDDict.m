
classdef nDDict

    properties
        data               % Storing the actual data (multi-dimensional matrix or cell array)
        axis = nDDictAxis  % 1xNdims - array of nDDictAxis classes for each axis. Ndims = ndims(data)
        meta = struct;     % Metadata about stuff that's stored in data
    end
    
    
    methods

        function obj = nDDict(data,axis,meta)
            if exist('data','var')
                obj.data = data;
            end
            
            if exist('axis','var')
                obj.axis = axis;
            end
            
            if exist('meta','var')
                obj.meta = meta;
            end
            
        end
        
        function [obj] = reset(obj)
            obj.data = [];
            obj.axis = nDDictAxis;
            obj.meta = struct;
        end

        function [obj2, ro] = subset(obj,varargin)
            % Define variables and check that all dimensions are consistent
            % ro - if regular expressions are used, returns the index
            % values discovered by the regular expression.
            
            % Verify that size of obj is correct
            checkDims(obj);
            
            % Validate that "selection" is the right size
            selection = varargin(:);
            Na = length(obj.axis);
            Nd = ndims(obj.data);
            
            % Fill out selection with empties if too short. For example, say
            % size(obj.data) is MxNx1x1. The user might enter a selection with
            % length(selection) = 2, such as selection = {[A],[B]}. In this
            % case, we convert it from to selection = {[A],[B],[],[]}.
            Ns = length(selection);
            if Ns < Na && Ns >= Nd
                selection2 = repmat({[]},1,Na);
                selection2(1:Ns) = selection;
                selection = selection2;
                clear selection2
            end
            
            % If Ns is still wrong dimensions, return error
            Ns = length(selection);
            if Ns ~= Na
                error('Number of inputs must match dimensionality of nDDict.data');
            end
            
            % Convert selection to index if using regular expressions
            ro = {};
            for i = 1:Ns
                if ischar(selection{i})
                    ro{i}(1,:) = obj.axis(i).values;
                    [selection{i} ro{i}(2,:)] = regex_lookup(obj.axis(i).values, selection{i});
                    
                end
            end
            
            % Initialize
            sz = size(obj);
            obj2 = obj;             % Create new class of same type as original
            obj2 = obj2.reset;
            obj2.meta = obj.meta;
            
            % First update each axis with the corresponding selection and
            % convert empty cells to code for full range.
            for i = 1:Ns

                if isempty(selection{i})
                    selection{i} = 1:sz(i);
                end
                
                obj2.axis(i) = obj.axis(i);       % Import axis information
                obj2.axis(i).values = obj.axis(i).values(selection{i});   % Overwrite values field; leave everything else the same.
            end
            
            % Update the data
            obj2.data = obj.data(selection{:});
            
            % Corrects number of axes. The above code automatically
            % converts obj.data from MxNx1x1 to MxN, whereas axis will stay
            % as MxNx1x1 (e.g. length of 4). Thus, fixAxes corrects this.
            % This should no longer be necessary!!!
            % obj2 = fixAxes(obj2);
            
        end
        
        function obj = importAxisNames(obj,ax_names)
            Nd = ndims(obj.data);
            if nargin < 2
                ax_names = cellfun(@num2str,num2cell(1:5),'UniformOutput',0);
            end
            
            if length(ax_names) ~= Nd
                error('Mismatch between number of axis names supplied and number of dimensions in dataset'); end

            for i = 1:ndims(obj.data)
                obj.axis(i).name = ax_names{i};
            end
        end
        
        function obj = importMeta(obj,meta_struct)
            obj.meta = meta_struct;
        end
        
        function [selection_out, startIndex] = findaxis(obj,str)
            % Returns the index of the axis with name matching str
            allnames = {obj.axis.name};
            [selection_out, startIndex] = regex_lookup(allnames, str);
        end
        
        function obj = mergeDims(obj,dims2pack)
            
            obj.checkDims;
            Nd2p = length(dims2pack);
            sz = size(obj.data);
            Nmerged = prod(sz(dims2pack));       % Final number entries in the merged dimension
            
            
            % % First, do axis names
            % Get cell array of all linearized axis values.
            inds = 1:Nmerged;
            [subs{1:Nd2p}] = ind2sub(sz(dims2pack),inds);
            temp = cell(1,Nd2p);
            for i = 1:Nd2p
                for j = 1:Nmerged
                    ax = obj.axis(dims2pack(i));
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
            obj.axis(dims2pack(1)).values = tempstr;
            obj.axis(dims2pack(1)).astruct.premerged_values = temp;
            
            
            % Give it a new axis name, reflecting the merger of all the
            % others
            allnames = {obj.axis(dims2pack).name};
            allnames = cat(1,allnames(:)',repmat({'_'},1,length(allnames)));
            obj.axis(dims2pack(1)).name = strcat(allnames{1:end-1});
            
            % Clear the remaining axes names
            for i = 2:Nd2p
                obj.axis(dims2pack(i)).name = ['Dim ' num2str(dims2pack(i))];
                obj.axis(dims2pack(i)).values = 1;
            end


            % % Now, work on obj.data
            dat = obj.data;
            obj.data = [];       % Clear obj.data for now, to save memory.
            
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
            
            % Undo the earlier permute, and put back into obj.data
            dat = ipermute(dat,[dims2pack,dims_remaining]);
            obj.data = dat;
        end
        
        function obj = packDim(obj,dim_src,dim_target)
            
            if nargin < 3; dim_target = dim_src; end
            checkDims(obj);
            
            % Make sure that obj.data is a cell array
            if ~iscell(obj.data); error('nDDict.data must be a cell array.'); end
            
            % Make sure that obj.data is a matrix
            temp = cellfun(@ismatrix,obj.data);
            if any(temp(:) ~= 1); error('nDDict.data must contain only matrices'); end
            % % To do: implement this so it works with cell arrays and nDDict
            % classes in the future too
            
            % Make sure target dimension in nDDict.data is a singleton
            temp = cellfun(@(x) size(x,dim_target),obj.data);
            if any(temp(:) ~= 1); error('Target dimension in nDDict.data needs to be size 1. Try reshaping contents of nDDict.data or choosing a different target dimension.'); end
            clear sz_targets
            
            % Bring chosen dimension to the front. Thus, we will be
            % merging along rows.
            Nd = ndims(obj.data);
            obj.data = permute(obj.data,[dim_src, 1:dim_src-1, dim_src+1:Nd]);
            
            % Temporarily linearize all other dimensions.
            sz0 = size(obj.data);
            obj.data = reshape(obj.data,sz0(1),prod(sz0(2:end)));
            
            % Add NaNs where needed
                % Note: to understand what this is doing, it really helps
                % to draw a diagram!
            sz = size(obj.data);
            empties = cellfun(@isempty,obj.data);    % 2D matrix with 1's marking empty cells
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
                        obj.data{curr_bad(i),bi(j)} = NaN*ones(size(obj.data{curr_good(1),bi(j)}));
                    end
                end
            end
            
            % Check that sizes and dimensionalities are compatible
            data_ndims = cellfun(@ndims,obj.data,'UniformOutput',true);
            if any(any(data_ndims ~= repmat(data_ndims(1,:),sz(1),1),1),2)
                error('Dimensions of nDDict.data not uniform along packing dimensions.');
            end
           
            
            data_sz = cellfun(@size,obj.data,'UniformOutput',false);
            data_sz_firsts = repmat(data_sz(1,:),sz(1),1);
            myfunc = @(x,y) any(x(:) ~= y(:));
            bool_size_mismatch = cellfun(myfunc,data_sz,data_sz_firsts);
            if any(bool_size_mismatch(:))
                error('Sizes of nDDict.data are not uniform along packing dimension. (This usually results form trying to combine populations with different numbers of cells.');
            end
            
            % Now, pack the data.
%             for j = 1:sz(2)
%                 % permute so that contents of each cell has the selected
%                 % dimenison along dimension 1
%                 for i = 1:sz(1)
%                     dat_curr = obj.data{i,j};
%                     Ndimdat = ndims(dat_curr);
%                     dat_curr = permute(dat_curr,[dim, 1:dim-1, dim+1:Ndimdat]); %whos dat_curr
%                     % dat_curr = permute(dat_curr,[2:dim, 1, dim+1:Ndimdat]);    %whos dat_curr       % Uncomment this stuff to make sure it permutes back
%                     obj.data{i,j} = dat_curr;
%                 end
%                 
%                 % Concatenate the data along dimension 1
%                 data_new{j} = cat(1,obj.data{j,:});
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
                obj.data{1,j} = cat(dim_target,obj.data{:,j});
                obj.data(2:end,j) = cell(sz(1)-1,1);     % Set remainder to empty
            end
            
            obj.data = obj.data(1,:);         % Keep only 1st dimension;
            sz(1) = 1;
            sz0(1) = 1;

            % Lastly, restore original dimensions
            % of nDDict.data
            obj.data = reshape(obj.data,sz0);
            obj.data = permute(obj.data,[2:dim_src, 1, dim_src+1:Nd]);
            
            % Also, update obj.axis
%             obj.axis(dim_src).values = 1;
            obj = setAxisDefaultNames(obj,dim_src);
            obj.axis(dim_src).name = ['Dim ' num2str(dim_src)];
            
        end
        
        function out = getaxisinfo(obj)
            % If no output arguments, prints axis info to the screen. If
            % output arguments are supplied, returns this information as a
            % strin.g
            
            if nargout > 0
                out = '';
            end
            
            for i = 1:length(obj.axis)
                
                out1 = obj.axis(i).getaxisinfo;
                temp = '';
                
                if nargout > 0
                    out = [out, temp, out1, '; ' ];
                else
                    temp = ['Axis ', num2str(i), ': '];
                    fprintf([temp, out1, '\n']);
                end
            end
            
            % Lastly output a summary of dimensionality comparing nDDict.axis
            % and nDDict.data. These should match up.
            if nargout == 0
                fprintf(['nDDict.axis dimensionality ' num2str(cellfun(@length,{obj.axis.values})) '\n']);
                fprintf(['nDDict.data dimensionality ' num2str(size(obj.data)) '\n']);
            end
        end
        
        % % % % % % % % % % % HOUSEKEEPING FUNCTIONS % % % % % % % % % % %
        function obj = fixAxes(obj)
            % This function forces the nDDict axis data to be updated to
            % match the dimensions of the data structure.
            % The convention of nDDict is to have always follow MATLAB
            % conventions for dimensionality. Thus, the size(obj.data)
            % command is used to determine the correct number of axes, and
            % axis is adjusted to match, adding or removing dimensions as
            % needed. If you are getting errors when running checkDims, you
            % should run this command.
            %
            % The one exception to MATLAB conventions is that there are
            % allowed to be more axes than there are dimensions in obj.data
            % as long as the number of entries in each of these axes is 1.

            Nd = ndims(obj.data);
            Na = length(obj.axis);

            % Sweep through all axes and make sure dimensions are correct.
            % Add new axes if needed, up to Nd.
            for i = 1:Nd
                obj = setAxisDefaultNames(obj,i);  % Sets axis #i to the default name
            end

            % Trim away excess values in axes
            if Na > Nd
                for i = Nd+1:Na
                    if length(obj.axis(i).values) > 1
                        obj.axis(i).values = obj.axis(i).values(1);
                        fprintf(['Extra values found in axis #' num2str(i) ' ' obj.axis(i).name '. Trimming \n']);
                    end
                end
                
            end
            
        end
        
        
        function checkDims(obj)
            % We enforce that size(obj.data) must always match up to 
            % length(obj.axis(i).values) for all i. We allow there to be
            % more axis than ndims(obj.data), but only if these axes have
            % lengths of 1.
            
            % Note, fixAxes fixes everything automatically.
            % Only call checkDims if you want to
            % be alerted to mismatches, but not to correct them. Use fixAxes to
            % automatically correct everything.
            
            sza = arrayfun(@(x) length(x.values),obj.axis);
            szd = size(obj.data);
            
            Nd = ndims(obj.data);
            Na = length(obj.axis);

            if Nd > Na
                error('checkDims: Error found! Number of dimensions in nDDict.data does not equal number of axes');
            end

            % For all dimensions in obj.data
            for i = 1:Nd
                if sza(i) ~= szd(i)
                    fprintf(['checkDims: Error found! obj.data dimension ',num2str(i), ...
                        ' is ', num2str(szd(i)) , ...
                        '. But corresponding axis \"', obj.axis(i).name , ...
                        '\" has ', num2str(sza(i)) , ' values. Dimension mismatch. \nTry running nDDict.getaxisinfo and then nDDict.fixAxes. \n' ]);
                    error(' ');
                end
            end
            
            % For additional axes beyond ndims(obj.data)
            ind = sza > 1;
            if any(ind(Nd+1:Na))
                ind2 = find(ind);
                ind2 = ind2(ind2 > Nd);
                fprintf(['checkDims: Error found! ndims(obj.data)=' num2str(Nd) ' but axis obj.axis(' num2str(ind2) ').values has ' num2str(sza(ind2)) ' entries.\n']);
                error(' ');
            end
            
            % All good if make it to here
        end

        
        % % % % % % % % % % % OVERLOADED FUNCTIONS % % % % % % % % % % %
        
        function varargout = size(obj,varargin)
            % Returns size of obj. This function is basically the same as
            % running size(obj.data) except we base it off of the dimensions
            % of obj.axis rather than obj.data. This has the effect of
            % returning 1's if length(obj.axis) > ndims(obj.data)
            
            checkDims(obj);
            
            [varargout{1:nargout}] = size(obj.data,varargin{:});
            
            % If function is called in the form sz = size(obj) OR size(obj),
            % return the length of each axis.
            if nargout <= 1 && nargin == 1
                Na = length(obj.axis);
                sz = zeros(1,Na);
                for i = 1:Na
                    sz(i) = length(obj.axis(i).values);
                end
                if nargout == 1; varargout{1} = sz; end
            end
        end
        
        function Nd = ndims(obj)
            checkDims(obj);
            Nd = length(obj.axis);
        end
        
        function obj = permute(obj,order)
            checkDims(obj);
            obj.data = permute(obj.data,order);
            obj.axis = obj.axis(order);
        end
        
        
        function obj = transpose(obj)
            checkDims(obj);
            Nd = ndims(obj.data);
            
            if Nd > 2; error('Can only transpose data with at most 2 dimensions');
            end
            
            obj.data = obj.data';
            obj.axis([1,2]) = obj.axis([2,1]);        % Axis should always be at least length=2.
        end
        
        function obj = squeeze(obj)
            % This is just like MATLAB's normal squeeze command. However,
            % there is one key difference:
            % Normally, if squeeze operates on a 1xN matrix, it will leave
            % it as 1xN. This function forces it to always return as Nx1
            
            checkDims(obj);
            
            % If data is bigger than a matrix, squeeze out dimensions that
            % are of size 1.
            sz = size(obj.data);
            if length(sz) > 2
                ind = sz~=1;
                obj.axis = obj.axis(ind);

                % Now squeeze obj.data
                obj.data = squeeze(obj.data);         % Normal squeeze command

%                 % Lastly, if the result is a row vector, force it to be a
%                 % column vector
%                 if isvector(obj.data) && ~iscolumn(obj.data)
%                     obj.data = obj.data';
%                 end
            else
                % Otherwise, if data is a matrix, remove all axis beyond
                % the first two. These should only be size 1 (e.g. "name"
                % axes anyways)
%                 szA = cellfun(@length,{obj.axis.values});
%                 ind = szA~=1;
%                 ind(1:2) = true;
                obj.axis = obj.axis(1:2);
            end
            
            % Make sure everything is good before returning.
            obj = obj.fixAxes;
            checkDims(obj);
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




function obj = setAxisDefaultNames(obj,dim)
    % Sets obj.axis(i) to default values
    
    % Get desired size of dataset
    sz_dim = size(obj.data,dim);
    
    % If axis doesn't already exist, create it. Otherwise, copy existing.
    if length(obj.axis) < dim
        ax_curr = nDDictAxis;
    else
        ax_curr = obj.axis(dim);
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
    
    obj.axis(dim) = ax_curr;
end


function [selection_out, startIndex] = regex_lookup(vals, selection)
    if ~ischar(vals{1}); error('nDDict.axis.values must be strings when using regular expressions');
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
