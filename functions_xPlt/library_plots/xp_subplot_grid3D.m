

function hsg = xp_subplot_grid3D (xp,transpose_on, display_mode)
	% xp must be 1D, 2D, or 3D
    
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
    
    if transpose_on && ismatrix(xp)
        xp = xp.transpose;
    end
    
    % Parameters
    %subplot_grid_options = {'no_zoom'};
    subplot_grid_options = {};
    
    sz = size(xp);
    
    if ndims(xp.data) <= 2
        N1 = sz(1);
        N2 = sz(2);
        
        
            if display_mode == 1 
                h0 = gcf; ha0 = gca;
                h = figure('visible','off');
            else
                figure;
            end
            
            hsg = subplot_grid(N1,N2,subplot_grid_options{:});
            c=0;
            for i = 1:N1
                for j = 1:N2
                    c=c+1;
                    hsg.set_gca(c);
                    xp.data{i,j}(); 
%                     xp2 = xp.subset(i,j);
%                     title(strrep(xp2.getaxisinfo,'_',' '));
                end
            end
            
            % Do labels for rows
            rowstr = setup_axis_labels(xp.axis(1));
            hsg.rowtitles(rowstr);
            
            % Do labels for columns
            colstr = setup_axis_labels(xp.axis(2));
            hsg.coltitles(colstr);
            
            if display_mode == 1
                
                cdata = print(h,'-RGBImage');
                close(h);

                % Restore original axes and display image
                figure(h0); axes(ha0);
                imshow(cdata);
                
            end
        
        
    elseif ndims(xp.data) == 3
        N1 = sz(1);
        N2 = sz(2);
        N3 = sz(3);
        
        for i = 1:N1
            figure;
            hsg(i) = subplot_grid(N2,N3,subplot_grid_options{:});
            hsg(i).figplace(N1,i);
            hsg(i).figtitle([xp.axis(1).name ': ' xp.axis(1).getvaluestring(i)]);
            c=0;
            for j = 1:N2
                for k = 1:N3
                    c=c+1;
                    hsg(i).set_gca(c);
                    xp.data{i,j,k}(); 
%                     xp2 = xp.subset(i,j);
%                     title(strrep(xp2.getaxisinfo,'_',' '));
                end
            end
            
            % Do labels for rows
            rowstr = setup_axis_labels(xp.axis(2));
            hsg(i).rowtitles(rowstr);
            
            % Do labels for columns
            colstr = setup_axis_labels(xp.axis(3));
            hsg(i).coltitles(colstr);
            
            
        end
        
    end
    
    
    
    
end

function outstr = setup_axis_labels(xpa)
    vals = xpa.getvaluescellstring;
    vals = strrep(vals,'_',' ');
    outstr = cell(size(vals));
    for j = 1:length(outstr)
        outstr{j} = {'',vals{j}};
    end
    outstr{round(end/2)}{1} = strrep(xpa.name,'_',' ');
end
