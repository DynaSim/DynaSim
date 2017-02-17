

function xp_handles_newfig (xp)
    % xp must be 1D
    
    for i = 1:length(xp.data)
        figure; xp.data{i}();
    end
end


