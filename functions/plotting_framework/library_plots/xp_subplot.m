

function xp_subplot (xp)
    % xp must be 1D
    
    N = length(xp.data);
    
    h = figure('Position',[440    52   617   746]);
    for i = 1:N
        
        figure(h);
        subplot(N,1,i); 
        xp.data{i}(); 
    end
end


