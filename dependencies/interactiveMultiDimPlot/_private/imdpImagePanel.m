function imdpImagePanel(hObject, eventdata, handles)

if hObject.Value % turn on
  hObject.String = 'Image [on]';
  createImagePanel()
else % turn off
  hObject.String = 'Image [off]';
  if isvalid(handles.ImagePanel.handle)
%     hObject.Units = 'normalized';
%     handles.ImagePanel.lastPos = hObject.Position;
    close(handles.ImagePanel.handle)
  else % figure was already closed, so reopen
    hObject.Value = 1;
    hObject.String = 'Image [on]';
    createImagePanel()
  end
end

% Update handles structure
guidata(hObject, handles);

  function createImagePanel()
    hFig = figure('Name','Image Panel','NumberTitle','off');
%     if isfield(handles.ImagePanel, 'lastPos')
%       hFig.Units = 'normalized';
%       hFig.Position = handles.ImagePanel.lastPos;
%     end
    axes(hFig, 'Position', [0 0 1 1]);
    handles.ImagePanel.handle = hFig;
  end

end