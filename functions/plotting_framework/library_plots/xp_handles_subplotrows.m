

function xp_handles_subplotrows (xp)
    % xp must be 1D
    
    N = length(xp.data);
    
    h = figure('Position',[440    52   617   746]);
    for i = 1:N
        cdata = xp.data{i}(); 
        figure(h);
        subplot(N,1,i); 
        imshow(cdata);
    end
end


