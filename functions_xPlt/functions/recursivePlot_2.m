
function varargout = recursivePlot_2(xp,function_handles,dimensions,function_handle_arguments)
%% varargout = recursivePlot(xp,function_handles,dimensions,function_handle_arguments)
%     Purpose: Takes in multidimensional data and a set of function handles for
%     plotting. Function handles operate on low-dimensional subspaces. 
%     recursivePlot assigns the high dimensional data in xp to the plotting
%     functions until all dimensions have been used up.
% 
%     Forms:
%         varargout = recursivePlot_2(xp,function_handles,dimensions)
%         varargout = recursivePlot_2(xp,function_handles,dimensions,function_handle_arguments)
% 
%     Inputs:
%       xp                          : xPlt structure with Nd dimensions
% 
%       function_handles            : Cell array of function handles. See below
%                                     for function handle specifications.
% 
%       dimensions                  : Cell array of vectors. Each vector corresponds
%                                     to a function handle. Specifies which dimensions
%                                     in xp to assign to each function handle.
% 
%       function_handle_arguments   : (optional) cell array of cell arrays, each containing
%                                     additional arguments to pass to each function
%                                     handle. There should be one cell array for
%                                     each function handle.
% 
%     Details:
%       The function_handles cell array should point to functions of the form:
%             varargout = function func (xp,varargin)
%       where xp is an xPlt object.
% 
% recursivePlot does the following
%     1. Works with the first 
% 
%     Examples:
%         See demos file.


%     function_handle_arguments - cell array of argument cell arrays to pass to
%     function_handles. Must have one cell array for each function_handle
%     passed. Use empty cell arrays for no arguments. E.g.
%     function_handle_arguments = { {}, {}, {1} } for nothing, nothing, and
%     1 as arguments.
    
    if nargin < 4
        function_handle_arguments = cell(size(function_handles));
        for i = 1:length(function_handle_arguments)
            function_handle_arguments{i} = {};
        end
    end

    sz = size(xp);
    if length(function_handles) > 1
        
        % xp_temp = permute(xp, [dimensions{1}, cell2mat(dimensions{2:end})]);
        
        % Set up a new xPlt object with dimensions matching current
        % function handle
        xp2 = xp.reset;                         % Create a new xPlt object        
        sz2 = sz(dimensions{1});
        ndims2 = length(dimensions{1});
        total_calls = prod(sz2);
        
        selection_curr = cell(1,length(sz));
        
        clear output_indices
        [xp2_indices{1:ndims2}] = ind2sub(sz2, 1:total_calls);      % xp2 will loop through chosen dimensions
        
        xp2_indices = num2cell(cat(1, xp2_indices{:})');
        
        xp_indices = cell(1, length(xp.axis));                      % Initialize xp_indices
        
        for call = 1:total_calls
            
            xp_indices(dimensions{1}) = xp2_indices(call, :);
            
            xp2.data{xp2_indices{call, :}} = @() recursivePlot_2(xp.subset(xp_indices{:}),...
                function_handles(2:end),dimensions(2:end),function_handle_arguments(2:end));
            
        end
        
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
            
        