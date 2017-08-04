function xp_scatter(xp)
    
    xp_size = sort(size(xp), 2, 'descend');
    if xp_size(1) ~= 1
        error('xp_scatter can only be used with a scalar xp object.')
    end
    
    xp_mat_size = sort(size(xp.data{1}), 2, 'descend');
    if xp_mat_size(2) ~= 2 || length(xp_mat_size > 1) > 2
        error('xp_scatter can only be used with a scalar xp object whose data is n x 2.')
    end

    scatter(xp.data{1}(:, 1), xp.data{1}(:, 2))
    
end