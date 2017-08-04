function obj = xp_matrix_transpose(obj)

obj.data = cellfun(@(x) x', obj.data, 'UniformOutput', 0);

meta = obj.meta;
meta.matrix_dim_1 = obj.meta.matrix_dim_2;
meta.matrix_dim_2 = obj.meta.matrix_dim_1;
obj = obj.importMeta(meta);

end