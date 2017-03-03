function varargout = imdpPlot(hObject, eventdata, handles)

nViewDims = handles.PlotPanel.nViewDims;
viewDims = handles.PlotPanel.viewDims;
% nAxDims = handles.PlotPanel.nAxDims;
hFig = handles.PlotPanel.figHandle;
hAx = handles.PlotPanel.axHandle;
mdData = handles.mdData;
dimNames = mdData.dimNames;

if ~isValidFigHandle(hFig)
  return
end

lFontSize = 14;
lMarkerSize = 16; %legend marker size

figure(hFig); % set hFig for gcf

if isfield(handles.data,'Label')
  colors = cat(1,handles.PlotPanel.Label.colors{:});
  markers = handles.PlotPanel.Label.markers;
  groups = handles.PlotPanel.Label.names;
  plotVarNum = handles.PlotPanel.Label.varNum;
  plotLabels = mdData.data{plotVarNum};
end

%% TODO: check axis order correct **********************************************
switch nViewDims
  case 1
    % 1 1d pane
    make1dPlot(hAx)
  case 2
    % 1 2d pane
    plotDims = find(viewDims);
    make2dPlot(hAx, plotDims)
  case 3
    % 3 2d panes + 1 3d pane = 4 subplots
    plotDims = find(viewDims);
    
    % 2d plots
    plotDims2d = combnk(plotDims,2);
    for iAx = 1:3
      ax2d = hAx(iAx);
      make2dPlot(ax2d, plotDims2d(iAx,:));
    end
    
    % 3d plot
    ax3d = hAx(iAx+1);
    if isgraphics(ax3d) && isempty(get(ax3d,'Children'))
      make3dPlot(ax3d, plotDims)
    end
  case 4
    % 6 2d panes + 4 3d pane = 10 subplots
  case 5
    % 10 2d panes + 10 3d pane = 20 subplots
  case 6
    % 15 2d panes = 15 subplots
  case 7
    % 21 2d panes = 21 subplots
  case 8
    % 28 2d panes = 28 subplots
end

% hFig.hide_empty_axes;

if nargout > 0
  varargout{1} = handles;
end

%% Sub functions
  function make3dPlot(hAx, plotDims)
    % x dim is plotDims(1)
    % y dim is plotDims(2)
    % z dim is plotDims(2)
    
    axes(hAx)
    
    sliceInd = handles.PlotPanel.axInd;
    sliceInd = num2cell(sliceInd);
    [sliceInd{plotDims}] = deal(':');
    
    % Get grid
    [x,y,z] = meshgrid(mdData.dimVals{plotDims(1)}, mdData.dimVals{plotDims(2)}, mdData.dimVals{plotDims(3)});
    g = plotLabels(sliceInd{:});
    
    % Linearize grid
    x = x(:)';
    y = y(:)';
    z = z(:)';
    g = g(:)';
    
    % Remove empty points
    emptyCells = cellfun(@isempty,g);
    x(emptyCells) = [];
    y(emptyCells) = [];
    z(emptyCells) = [];
    g(emptyCells) = [];
    
    plotData.x = x;
    plotData.xlabel = dimNames{plotDims(1)};
    
    plotData.y = y;
    plotData.ylabel = dimNames{plotDims(2)};
    
    plotData.z = z;
    plotData.zlabel = dimNames{plotDims(3)};
    
    plotData.g = g;
    
    plotData.clr = [];
    plotData.sym = '';
    for grp = unique(plotData.g)
      gInd = strcmp(groups, grp);
      thisClr = colors(gInd,:);
      thisSym = markers{gInd};
      plotData.clr(end+1,:) = thisClr;
      plotData.sym = [plotData.sym thisSym];
    end
    
    % Marker Size
    if handles.autoSizeMarkerCheckbox.Value %auto size marker
      set(gca,'unit', 'pixels');
      pos = get(gca,'position');
      axSize = pos(3:4);
      markerSize = min(axSize) / max([length(plotData.x), length(plotData.y), length(plotData.z)]);
      set(gca,'unit', 'normalized');
      plotData.siz = markerSize;
    else %manual size marker
      markerSize = handles.markerSizeSlider.Value;
      plotData.siz = markerSize;
    end
    
    % Set MarkerSize Slider Val
    if isfield(handles.PlotPanel, 'sliderH')
      handles.PlotPanel.sliderH.Value = markerSize;
      imdpMarkerSizeSliderCallback(handles.PlotPanel.sliderH,[])
    end
    
    scatter3dPlot(plotData);
    
%     % Rescale ylim
%     ylims = get(hAx,'ylim');
%     set(hAx, 'ylim', [ylims(1)- 0.05*range(ylims) ylims(2)+0.05*range(ylims)]);
  end
  
  function make2dPlot(hAx, plotDims)
    % x dim is plotDims(1)
    % y dim is plotDims(2)
    
    axes(hAx)
    
    sliceInd = handles.PlotPanel.axInd;
    sliceInd = num2cell(sliceInd);
    [sliceInd{plotDims}] = deal(':');
    
    % Get grid
    [x,y] = meshgrid(mdData.dimVals{plotDims(1)}, mdData.dimVals{plotDims(2)});
    g = plotLabels(sliceInd{:});
    
    % Linearize grid
    x = x(:)';
    y = y(:)';
    g = g(:)';
    
    % Remove empty points
    emptyCells = cellfun(@isempty,g);
    x(emptyCells) = [];
    y(emptyCells) = [];
    g(emptyCells) = [];
    
    plotData.x = x;
    plotData.xlabel = dimNames{plotDims(1)};
    
    plotData.y = y;
    plotData.ylabel = dimNames{plotDims(2)};
    
    plotData.g = g;
    
    plotData.clr = [];
    plotData.sym = '';
    for grp = unique(plotData.g)
      gInd = strcmp(groups, grp);
      thisClr = colors(gInd,:);
      thisSym = markers{gInd};
      plotData.clr(end+1,:) = thisClr;
      plotData.sym = [plotData.sym thisSym];
    end
    
    % Marker Size
    if handles.autoSizeMarkerCheckbox.Value %auto size marker
      set(gca,'unit', 'pixels');
      pos = get(gca,'position');
      axSize = pos(3:4);
      markerSize = min(axSize) / max(length(plotData.x), length(plotData.y));
      set(gca,'unit', 'normalized');
      plotData.siz = markerSize;
    else %manual size marker
      markerSize = handles.markerSizeSlider.Value;
      plotData.siz = markerSize;
    end
    
    % Set MarkerSize Slider Val
    if isfield(handles.PlotPanel, 'sliderH')
      handles.PlotPanel.sliderH.Value = markerSize;
      imdpMarkerSizeSliderCallback(handles.PlotPanel.sliderH,[])
    end
    
    scatter2dPlot(plotData);
    
    % Rescale ylim
    try
      ylims = get(hAx,'ylim');
      set(hAx, 'ylim', [ylims(1)- 0.05*range(ylims) ylims(2)+0.05*range(ylims)]);
    end
  end

  function make1dPlot(hAx)
    axes(hAx)
    plotDim = find(viewDims);
    sliceInd = handles.PlotPanel.axInd;
    sliceInd = num2cell(sliceInd);
    sliceInd{plotDim} = ':';
    
    plotData.xlabel = dimNames{plotDim};
    plotData.x = mdData.dimVals{plotDim};
    plotData.y = zeros(length(plotData.x),1);
    plotData.ylabel = '';
    plotData.g = plotLabels(sliceInd{:});
    plotData.g = plotData.g(:)';
    
    % Remove empty points
    emptyCells = cellfun(@isempty,plotData.g);
    plotData.x(emptyCells) = [];
    plotData.y(emptyCells) = [];
    plotData.g(emptyCells) = [];
    
    plotData.clr = [];
    plotData.sym = '';
    for grp = unique(plotData.g)
      gInd = strcmp(groups, grp);
      thisClr = colors(gInd,:);
      thisSym = markers{gInd};
      plotData.clr(end+1,:) = thisClr;
      plotData.sym = [plotData.sym thisSym];
    end
    
    % Marker Size
    if handles.autoSizeMarkerCheckbox.Value %auto size marker
      set(gca,'unit', 'pixels');
      pos = get(gca,'position');
      axSize = pos(3:4);
      markerSize = min(axSize) / length(plotData.x);
      set(gca,'unit', 'normalized');
      plotData.siz = markerSize;
    else %manual size marker
      markerSize = handles.markerSizeSlider.Value;
      plotData.siz = markerSize;
    end
    
%     % Set MarkerSize Slider Val
%     if isfield(handles.PlotPanel, 'sliderH')
%       handles.PlotPanel.sliderH.Value = markerSize;
%       imdpMarkerSizeSliderCallback(handles.PlotPanel.sliderH,[])
%     end
    
    scatter2dPlot(plotData);

    set(gca,'YTick', []);
  end

  function scatter2dPlot(plotData)
    try
      gscatter(plotData.x,plotData.y,categorical(plotData.g),plotData.clr,plotData.sym,plotData.siz,'off',plotData.xlabel,plotData.ylabel)
      if handles.MainPanel.legendBool
        uG = unique(plotData.g);
        [lH,icons] = legend(uG); % TODO: hide legend before making changes   

        % Increase legend width
    %     lPos = lH.Position;
    %     lPos(3) = lPos(3) * 1.05; % increase width of legend
    %     lH.Position = lPos;

        [icons(1:length(uG)).FontSize] = deal(lFontSize);
        [icons(1:length(uG)).FontUnits] = deal('normalized');

        shrinkText2Fit(icons(1:length(uG)))

        [icons(length(uG)+2:2:end).MarkerSize] = deal(lMarkerSize);
        
%         legend(gca,'boxoff')
%         legend(gca,'Location','SouthEast')
      end
    end
  end

  function scatter3dPlot(plotData)
    %     [uniqueGroups, uga, ugc] = unique(group);
    %     colors = colormap;
    %     markersize = 20;
    %     scatter3(x(:), y(:), z(:), markersize, colors(ugc,:));
    
    try
      [~, ~, groupInd4color] = unique(plotData.g);
      
%       plotData.sym

      scatter3(plotData.x, plotData.y, plotData.z, plotData.siz, plotData.clr(groupInd4color,:), '*');
      
      xlabel(plotData.xlabel)
      ylabel(plotData.ylabel)
      zlabel(plotData.zlabel)
      
%       if handles.MainPanel.legendBool
%         uG = unique(plotData.g);
%         [lH,icons] = legend(uG); % TODO: hide legend before making changes   
% 
%         % Increase legend width
%     %     lPos = lH.Position;
%     %     lPos(3) = lPos(3) * 1.05; % increase width of legend
%     %     lH.Position = lPos;
% 
%         [icons(1:length(uG)).FontSize] = deal(lFontSize);
%         [icons(1:length(uG)).FontUnits] = deal('normalized');
% 
%         shrinkText2Fit(icons(1:length(uG)))
% 
%         [icons(length(uG)+2:2:end).MarkerSize] = deal(lMarkerSize);
%         
% %         legend(gca,'boxoff')
% %         legend(gca,'Location','SouthEast')
%       end
    end
  end

  function shrinkText2Fit(txtH)
    for iTxt=1:length(txtH)
      % Check width
      ex = txtH(iTxt).Extent;
      bigBool = ( (ex(1) + ex(3)) > 1 );
      while bigBool
        txtH(iTxt).FontSize = txtH(iTxt).FontSize * 0.99;
        ex = txtH(iTxt).Extent;
        bigBool = ( (ex(1) + ex(3)) > 1 );
      end
    end
  end

end