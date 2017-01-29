

function cdata = xp_images_subplotcols (xp)
    % xp must be 1D
    
    N = length(xp.data);
    h = figure('visible','off','Position',[ 440   659   497   139]);
    %h = figure('Position',[ 440   659   497   139]);
    for i = 1:N
        subplot(1,N,i); xp.data{i}();
    end
    %pause(1);
    cdata = print(h,'-RGBImage');
    %print(h,['fig_' num2str(i) '_' num2str(randn) '.png'],'-dpng')
    %pause(1);
    close(h);
end


