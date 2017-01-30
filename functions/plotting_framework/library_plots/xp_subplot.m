

function xp_subplot (xp,transpose_on)
    % xp must be 1D or 2D
    
    if nargin < 2
        transpose_on = 0;
    end
    
    if transpose_on
        xp = xp.transpose;
    end
    
    sz = size(xp);
    N1 = sz(1);
    N2 = sz(2);
    
    h = figure('Position',[440    52   617   746]);
    c=0;
    for i = 1:N1
        for j = 1:N2
            c=c+1;
            figure(h);
            subplot(N1,N2,c); 
            xp.data{i,j}(); 
%             xp2 = xp.select(i,j);
%             title(xp2.getaxisinfo);
        end
    end
end


