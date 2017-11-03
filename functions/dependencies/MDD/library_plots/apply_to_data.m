function xp_out = apply_to_data(xp, function_handle, varargin)

xp_out = feval(function_handle, xp.data{:}, varargin{:});

end