function obj_out = mean_over_axis(obj, axis_dim, function_handle, varargin)

    if nargin < 3, function_handle = []; end
    
    if isempty(function_handle), function_handle = @nanmean; end
    
    if nargin < 4, varargin = {}; end
    
    if ischar(axis_dim)
        dim_string = axis_dim;
        axis_dim = obj.findaxis(dim_string);
        if ~isscalar(axis_dim) || isempty(axis_dim)
            error('Multiple or zero dimensions matching %s.', dim_string)
        end
    end
    
    if ~isscalar(axis_dim) || isempty(axis_dim) || axis_dim == 0
        error('Dimension to normalize must be a nonempty, nonzero scalar.')
    end
   
    %% Packing axis to take mean over.
    
    mean_dim = obj.lastNonSingletonDim + 1;
    
    obj_out = obj.packDim(axis_dim, mean_dim);
    
    %% Taking mean.
    
    varargin{end + 1} = mean_dim;
    
    obj_out.data = cellfun(@(x) feval(function_handle, x, varargin{:}), obj_out.data, 'UniformOutput', 0);
    
    obj_out.meta = rmfield(obj_out.meta, ['matrix_dim_', num2str(mean_dim)]);
    
    obj_out = obj_out.permute([1:(axis_dim - 1) (axis_dim + 1):ndims(obj_out) axis_dim]);
    
    obj_out.axis(end) = [];
    
end