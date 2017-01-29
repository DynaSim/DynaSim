
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
            
            % Initialize
            sz = size(xp);
            xp2 = xPlt;
            xp2.meta = xp.meta;
            
            % First update each axis with the corresponding selection
            for i = 1:length(selection)
                if all(cellfun(@length,selection(i:end)) == 1)
                    continue
                    % Stop adding axes if all remainding dims are 1
                    % Matlab automatically squeezes all trailing
                    % dimensions, so we will follow suit.
                end
                if isempty(selection{i})
                    selection{i} = 1:sz(i);
                end
                
                xp2.axis(i) = xp.axis(i);       % Import axis information
                xp2.axis(i).values = xp.axis(i).values(selection{i});   % Overwrite values field; leave everything else the same.
            end

            % Lastly update the data
            xp2.data = xp.data(selection{:});
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
            error('incomplete');
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
        
        
        
        function getaxisinfo(xp)
            for i = 1:length(xp.axis)
                fprintf(['Axis ', num2str(i), ': ']);
                xp.axis(i).getaxisinfo;
            end
        end
        
        % % % % % % % % % % % OVERLOADED FUNCTIONS % % % % % % % % % % %
        
        function sz = size(xp)
            % Overrides normal size command.
            for j = 1:length(xp.axis)
                sz(j) = length(xp.axis(j).values);
            end
        end
        
        
        function xp = permute(xp,order)
            xp.data = permute(xp.data,order);
            xp.axis = xp.axis(order);
        end
        
        function xp = squeeze(xp)
            % This is just like MATLAB's normal squeeze command. However,
            % there is one key difference:
            % Normally, if squeeze operates on a 1xN matrix, it will leave
            % it as 1xN. This function forces it to always return as Nx1
            
            checkDims(xp);
            
            % Remove axes that have dimensionality of 1
            sz = size(xp.data);
            xp.axis = xp.axis(sz~=1);
            
            % Now squeeze xp.data
            xp.data = squeeze(xp.data);         % Normal squeeze command
            
            % Lastly, if the result is a row vector, force it to be a
            % column vector
            if isvector(xp.data) && ~iscolumn(xp.data)
                xp.data = xp.data';
            end
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




function checkDims(xp)
    sz = size(xp.data);
    Nd = ndims(xp.data);
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
                '\" has ', num2str(Nvalues_in_axis) , ' elements. Dimension mismatch. \nTry running xPlt.getaxisinfo. \n' ]);
            error(' ');
        end
    end
end
