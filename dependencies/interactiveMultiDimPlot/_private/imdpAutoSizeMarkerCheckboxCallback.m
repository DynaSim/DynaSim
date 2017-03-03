function imdpAutoSizeMarkerCheckboxCallback(hObject, eventdata, handles)

% Enable/Disable marker size slider based on state of autoSize checkbox
if hObject.Value % auto on
  handles.markerSizeSlider.Enable = 'off';
else % auto off
  handles.markerSizeSlider.Enable = 'on';
end

imdpPlot(hObject, eventdata, handles);

end