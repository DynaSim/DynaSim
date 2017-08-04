

function hxp = xp_subplot_grid3D (xp, op)
	% This handles 1D, 2D, or 3D xp data. 3D data is tiled across the
	% screen in different figures.
    
    hxp = struct; 
    
    if nargin < 2
        op = struct;
    end
    
    if isempty(op); op = struct; end;
    
    op = struct_addDef(op,'transpose_on',0);
    op = struct_addDef(op,'display_mode',0);
    op = struct_addDef(op,'subplotzoom_enabled',1);
    op = struct_addDef(op,'legend1',[]);
    op = struct_addDef(op,'max_legend',20);
            % Display_mode: 0-Just plot directly
                          % 1-Plot as an image (cdata)
                          % 2-Save to a figure file 
                          
    transpose_on = op.transpose_on;
    display_mode = op.display_mode;
    subplotzoom_enabled = op.subplotzoom_enabled;
    legend1 = op.legend1;
                          
    if verLessThan('matlab','8.4') && display_mode == 1; warning('Display_mode==1 might not work with earlier versions of MATLAB.'); end
    if transpose_on && ismatrix(xp)
        xp = xp.transpose;
    elseif transpose_on && ~ismatrix(xp.data)
        error('xp must be a matrix (e.g. ndims < 3) in order to use transpose');
    end
    
    % Parameters
    %subplot_grid_options = {'no_zoom'};
    subplot_grid_options = {};
    
    sz = size(xp);
    
    if ndims(xp.data) <= 2
        hxp = xp_subplot_grid (xp, display_mode, transpose_on);
        
        
    elseif ndims(xp.data) == 3
        N1 = sz(1);
        N2 = sz(2);
        N3 = sz(3);
        
        for i = 1:N1
            figure;
            hxp.hcurr(i) = subplot_grid(N2,N3,subplot_grid_options{:});
            if subplotzoom_enabled
                hxp.hcurr(i) = subplot_grid(N2,N3,subplot_grid_options{:});
            else
                hxp.hcurr(i) = subplot_grid(N2,N3,'no_zoom',subplot_grid_options{:});
            end
            
            if ~verLessThan('matlab','8.4'); hxp.hcurr(i).figplace(N1,i); end
            mytitle = [figformat_str(xp.axis(1).name) ': ' figformat_str(xp.axis(1).getvalues_cellstr{i})];
            hxp.hcurr(i).figtitle(mytitle);
            c=0;
            for j = 1:N2
                for k = 1:N3
                    c=c+1;
                    hxp.hcurr(i).set_gca(c);
                    hxp.hsub = xp.data{i,j,k}(); 
                    if j == 1 && k == 1 && ~isempty(legend1)
                        % Place a legend in the 1st subplot
                        legend(legend1{1:min(end,op.max_legend)});
                    end
                end
            end
            
            % Do labels for rows
            if ~strcmp(xp.axis(2).name(1:3),'Dim')          % Only display if its not an empty axis
                rowstr = setup_axis_labels(xp.axis(2));
                hxp.hcurr(i).rowtitles(rowstr);
            end
            
            % Do labels for columns
            if ~strcmp(xp.axis(3).name(1:3),'Dim')          % Only display if its not an empty axis
                colstr = setup_axis_labels(xp.axis(3));
                hxp.hcurr(i).coltitles(colstr);
            end
            
        end
        
    end
    
end

function outstr = setup_axis_labels(xpa)
    vals = xpa.getvalues_cellstr;
    vals = strrep(vals,'_',' ');
    outstr = cell(size(vals));
    for j = 1:length(outstr)
        outstr{j} = {'',vals{j}};
    end
    outstr{round(end/2)}{1} = strrep(xpa.name,'_',' ');
end
