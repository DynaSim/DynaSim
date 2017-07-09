
function obj = xp_abs(obj)

checkDims(obj);

obj.data = cellfun(@(x) abs(x), obj.data_pr, 'UniformOutput', 0);

end