
function varargout = recursivePlot_2(xp,function_handles,dimensions,function_handle_arguments)
    % function_handle_arguments - cell array of argument cell arrays to pass to
    % function_handles. Must have one cell array for each function_handle
    % passed. Use empty cell arrays for no arguments. E.g.
    % function_handle_arguments = { {}, {}, {1} } for nothing, nothing, and
    % 1 as arguments.
    
    if nargin < 4
        function_handle_arguments = cell(size(function_handles));
        for i = 1:length(function_handle_arguments)
            function_handle_arguments{i} = {};
        end
    end

    sz = size(xp);
    if length(function_handles) > 1
        
        % xp_temp = permute(xp, [dimensions{1}, cell2mat(dimensions{2:end})]);
        
        xp2 = xPlt;
        selection_curr = cell(1,length(sz));
        % Just cardcode in the various cases for dimensionality. It's
        % unlikely they will ever want to go above showing 4D in a single
        % plot.
        
        temp = size(xp);
        
        current_size = temp(dimensions{1});
        
        current_no_dims = length(dimensions{1});
        
        total_calls = prod(current_size);
        
        [output_indices{1:current_no_dims}] = ind2sub(current_size, 1:total_calls);
        
        output_indices = num2cell(cat(1, output_indices{:})');
         
        input_indices = cell(1, length(xp.axis));
        
        for call = 1:total_calls
            
            input_indices(dimensions{1}) = output_indices(call, :);
            
            xp2.data{output_indices{call, :}} = @() recursivePlot(xp.subset(input_indices{:}),...
                function_handles(2:end),dimensions(2:end),function_handle_arguments(2:end));
            
        end
        
%         switch length(dimensions{1})        
%             case 1                          % 1D
%                 dim1 = dimensions{1}(1);
%                 for i = 1:sz(dim1)
%                     selection_curr{dim1} = i;
%                         % Note: need to make it a row vector!
%                     xp2.data{i,1} = @() recursivePlot(xp.subset(selection_curr{:}),function_handles(2:end),dimensions(2:end),function_handle_arguments(2:end));
%                 end
%                 
%             case 2                          % 2D
%                 dim1 = dimensions{1}(1);
%                 dim2 = dimensions{1}(2);
%                 for i = 1:sz(dim1)
%                     for j = 1:sz(dim2)
%                         selection_curr{dim1} = i;
%                         selection_curr{dim2} = j;
%                         xp2.data{i,j} = @() recursivePlot(xp.subset(selection_curr{:}),function_handles(2:end),dimensions(2:end),function_handle_arguments(2:end));
%                     end
%                 end
%                 
%             case 3                          % 3D
%                 dim1 = dimensions{1}(1);
%                 dim2 = dimensions{1}(2);
%                 dim3 = dimensions{1}(3);
%                 for i = 1:sz(dim1)
%                     for j = 1:sz(dim2)
%                         for k = 1:sz(dim3)
%                             selection_curr{dim1} = i;
%                             selection_curr{dim2} = j;
%                             selection_curr{dim3} = k;
%                             xp2.data{i,j,k} = @() recursivePlot(xp.subset(selection_curr{:}),function_handles(2:end),dimensions(2:end),function_handle_arguments(2:end));
%                         end
%                     end
%                 end
%                 
%             case 4                          % 4D
%                 dim1 = dimensions{1}(1);
%                 dim2 = dimensions{1}(2);
%                 dim3 = dimensions{1}(3);
%                 dim4 = dimensions{1}(4);
%                 for i = 1:sz(dim1)
%                     for j = 1:sz(dim2)
%                         for k = 1:sz(dim3)
%                             for l = 1:sz(dim4)
%                                 selection_curr{dim1} = i;
%                                 selection_curr{dim2} = j;
%                                 selection_curr{dim3} = k;
%                                 selection_curr{dim4} = l;
%                                 xp2.data{i,j,k,l} = @() recursivePlot(xp.subset(selection_curr{:}),function_handles(2:end),dimensions(2:end),function_handle_arguments(2:end));
%                             end
%                         end
%                     end
%                 end
%                 
%                 
%         end
        
        % Update axes
        for i = 1:length(dimensions{1})
            dim_curr = dimensions{1}(i);
            xp2.axis(i).values = xp.axis(dim_curr).values;
            xp2.axis(i).name = xp.axis(dim_curr).name;
        end
        xp2 = xp2.fixAxes;
        
        [varargout{1:nargout}] = function_handles{1}(xp2,function_handle_arguments{1}{:});
        
    else
        [varargout{1:nargout}] = function_handles{1}(xp,function_handle_arguments{1}{:});
    end
end
            
        