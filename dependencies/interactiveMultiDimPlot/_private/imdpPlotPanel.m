function varargout = imdpPlotPanel(hObject, eventdata, handles)

nViewDims = handles.PlotPanel.nViewDims;
nViewDimsLast = handles.PlotPanel.nViewDimsLast;

% if ~isfield(handles.PlotPanel, 'figHandle') || ~isvalid(handles.PlotPanel.figHandle) || nViewDims ~= nViewDimsLast
if ~isValidFigHandle('handles.PlotPanel.figHandle') || nViewDims ~= nViewDimsLast
  
  % Make New Panel
%   if ~isfield(handles.PlotPanel, 'figHandle') || ~isvalid(handles.PlotPanel.figHandle)% && nViewDims == nViewDimsLast
  if ~isValidFigHandle('handles.PlotPanel.figHandle') % && nViewDims == nViewDimsLast
    hFig = figure('Name','Plot Panel','NumberTitle','off');
    handles.PlotPanel.figHandle = hFig;
    
    % Update handles structure
    guidata(hObject, handles);
    
    newPanelBool = true;
  else
    hFig = handles.PlotPanel.figHandle;
    newPanelBool = false;
  end
  
  % Add user data to figure
  hFig.UserData.mdData = handles.mdData;
  
  % Data cursor
  dm = datacursormode(hFig);
  dm.UpdateFcn = @imdpDataCursorCallback;
  
  %   if isfield(handles.PlotPanel, 'handle') && isvalid(handles.PlotPanel.figHandle)
  %     hFig = handles.PlotPanel.figHandle;
  %
  %     if ~ishandle(hFig)
  %       close(hFig.hfig)
  %       hFig = figure;
  %       axes(hFig)
  %     end
  %   else
  %     hFig = figure;
  %     axes(hFig)
  %     handles.PlotPanel.figHandle = hFig;
  %
  %     % Update handles structure
  %     guidata(hObject, handles);
  %   end
  
  % Update Panel
  clf(hFig) %clear fig
  gap = 0.1;
  marg_h = 0.1;
  marg_w = 0.1;
  switch nViewDims
    case 1
      % 1 1d pane
      %         axes(hFig)
      %       hspg = subplot_grid(1,'no_zoom', 'parent',hFig);
      hAx = tight_subplot2(1, 1, gap, marg_h, marg_w, hFig);
    case 2
      % 1 2d pane
      %         axes(hFig)
      %       hspg = subplot_grid(1,'no_zoom', 'parent',hFig);
      hAx = tight_subplot2(1, 1, gap, marg_h, marg_w, hFig);
    case 3
      % 3 2d panes + 1 3d pane = 4 subplots
      %       hspg = subplot_grid(2,2, 'parent',hFig);
      hAx = tight_subplot2(2, 2, gap, marg_h, marg_w, hFig);
    case 4
      % 6 2d panes + 4 3d pane = 10 subplots
      %       hspg = subplot_grid(2,5, 'parent',hFig);
      hAx = tight_subplot2(2, 5, gap, marg_h, marg_w, hFig);
    case 5
      % 10 2d panes + 10 3d pane = 20 subplots
      %       hspg = subplot_grid(3,7, 'parent',hFig); % 1 empty
      hAx = tight_subplot2(3, 7, gap, marg_h, marg_w, hFig);
    case 6
      % 15 2d panes = 15 subplots
      %       hspg = subplot_grid(3,5, 'parent',hFig);
      hAx = tight_subplot2(3, 5, gap, marg_h, marg_w, hFig);
    case 7
      % 21 2d panes = 21 subplots
      %       hspg = subplot_grid(3,7, 'parent',hFig);
      hAx = tight_subplot2(3, 7, gap, marg_h, marg_w, hFig);
    case 8
      % 28 2d panes = 28 subplots
      %       hspg = subplot_grid(4,7, 'parent',hFig);
      hAx = tight_subplot2(4, 7, gap, marg_h, marg_w, hFig);
    otherwise
      if newPanelBool
        fprintf('\nSelect at least 1 ViewDim to plot.\n')
      end
  end
  
  %   if exist('hspg', 'var')
  %     handles.PlotPanel.figHandle = hspg;
  %
  %     % Update handles structure
  %     guidata(hObject, handles);
  %
  %     % Plot
  %     imdpPlot(hObject, eventdata, handles);
  %   end
  
  if nViewDims > 0
    % Axis handle
    handles.PlotPanel.axHandle = hAx;
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Plot
    handles = imdpPlot(hObject, eventdata, handles);
  end
end

if nargout > 0
  varargout{1} = handles;
end

end
