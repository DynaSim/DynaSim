

function h=xp_subplot (xp,transpose_on,display_mode)
    % xp must be 1D or 2D
    
    if nargin < 2
        transpose_on = [];
    end
    
    if nargin < 3
        display_mode = [];
    end
    
    if isempty(transpose_on); transpose_on = 0; end
    if isempty(display_mode); display_mode = 0; end
            % Display_mode: 0-Just plot directly
                          % 1-Plot as an image (cdata)
                          % 2-Save to a figure file
    
    if transpose_on
        xp = xp.transpose;
    end
    
    
    
    if display_mode == 0
        
        h=xp_subplot1_noimage(xp);
        
    elseif display_mode == 1
        
        h=xp_subplot2_asimage(xp);

    end
end



function h=xp_subplot1_noimage(xp)
        
        sz = size(xp);
        N1 = sz(1);
        N2 = sz(2);
        
        % Creates a new figure with a bunch of subplots
        h = figure('Units','normalized','Position',[0,0,1,1]);
        c=0;
        for i = 1:N1
            for j = 1:N2
                c=c+1;
                figure(h);
                subplot(N1,N2,c); 
                xp.data{i,j}(); 
                xp2 = xp.subset(i,j);
                title(strrep(xp2.getaxisinfo,'_',' '));
            end
        end
end

function h0 = xp_subplot2_asimage(xp)

    sz = size(xp);
    N1 = sz(1);
    N2 = sz(2);

    % Creates an existing figure as an image
    h0 = gcf; ha0 = gca;
    %h = figure('visible','off','Position',[ 440   659   497   139]);
    h = figure('visible','off');
    c=0;
    for i = 1:N1
        for j = 1:N2
            c=c+1;
            subplot(N1,N2,c); 
            xp.data{i,j}(); 
            xp2 = xp.subset(i,j);
            title(strrep(xp2.getaxisinfo,'_',' '));
        end
    end
    cdata = print(h,'-RGBImage');
    close(h);

    % Restore original axes and display image
    figure(h0); axes(ha0);
    imshow(cdata);
end