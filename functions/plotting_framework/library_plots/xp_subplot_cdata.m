

function cdata = xp_subplot_cdata (xp)
    % xp must be 1D
    
    N = length(xp.data);
    h0 = gcf; ha0 = gca;
    h = figure('visible','off','Position',[ 440   659   497   139]);
    for i = 1:N
        subplot(1,N,i); xp.data{i}();
    end
    cdata = print(h,'-RGBImage');
    close(h);
    
    % Restore original axes
    figure(h0); axes(ha0);
    imshow(cdata);
    
end


