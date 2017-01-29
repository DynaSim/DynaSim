

function xp_handles_subplotrows (xp)
    % xp must be 1D
    
    N = length(xp.data);
    for i = 1:N
        subplot(N,1,i); xp.data{i}();
    end
end


