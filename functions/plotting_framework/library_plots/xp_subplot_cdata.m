

function cdata = xp_subplot_cdata (xp,transpose_on)
    % xp must be 1D
    
    if nargin < 2
        transpose_on = 0;
    end
    
    if transpose_on
        xp = xp.transpose;
    end
    
    sz = size(xp);
    N1 = sz(1);
    N2 = sz(2);
    
    N = length(xp.data);
    h0 = gcf; ha0 = gca;
    h = figure('visible','off','Position',[ 440   659   497   139]);
    c=0;
    for i = 1:N1
        for j = 1:N2
            c=c+1;
            subplot(N1,N2,c); 
            xp.data{i,j}(); 
%             xp2 = xp.subset(i,j);
%             title(xp2.getaxisinfo);
        end
    end
    cdata = print(h,'-RGBImage');
    close(h);
    
    % Restore original axes and display image
    figure(h0); axes(ha0);
    imshow(cdata);
    
end


