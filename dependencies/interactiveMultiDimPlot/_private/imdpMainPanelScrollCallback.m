function imdpMainPanelScrollCallback(figH, scrollData)

mousePos = get(figH, 'CurrentPoint');

handles = imdpHandlesFromFig(figH);

nVarSliders = length(handles.MainPanel.Handles.sH);

onSlider = 0;
for iSlider = 1:nVarSliders + 1
  if iSlider <= nVarSliders % var slider
    thisSliderPos = handles.MainPanel.Handles.sH{iSlider}.Position;
  else % marker size slider
    thisSliderPos = handles.markerSizeSlider.Position;
  end
  
  x = thisSliderPos(1);
  y = thisSliderPos(2);
  w = thisSliderPos(3);
  h = thisSliderPos(4);
  
  % corners of slider area clocwise from LR corner
  sliderX = [x,x,x+w,x+w];
  sliderY = [y,y+h,y+h,y];
  inBool = inpolygon(mousePos(1),mousePos(2),sliderX,sliderY);
  if inBool
    onSlider = iSlider;
  end
end

% fprintf('Scrolling over slider #: %i\n', onSlider)

if onSlider
  eventdata = scrollData.VerticalScrollAmount*scrollData.VerticalScrollCount;
  
  if onSlider <= nVarSliders % var slider
    hObject = handles.MainPanel.Handles.sH{onSlider};
    imdpSliderChangeCallback(hObject, eventdata, handles)
  else% marker size slider
    hObject = handles.markerSizeSlider;
    imdpMarkerSizeSliderCallback(hObject, eventdata, handles)
  end
end

end
