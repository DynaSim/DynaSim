

function xp_matrix_imagesc (xp)
    % xp must be 1x1 (e.g. zero dimensional)
    
    imagesc(xp.data{1}');
    %colorbar
    
end


