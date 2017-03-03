function imdpMarkerSizeSliderCallback(hObject, eventdata, handles)

if isnumeric(eventdata) %from scroll
  newVal = hObject.Value + eventdata * (hObject.SliderStep(2))*(hObject.Max-hObject.Min);
  if newVal > hObject.Max
    newVal = hObject.Max;
  elseif newVal < hObject.Min
    newVal = hObject.Min;
  end
  hObject.Value = newVal;
end

imdpPlot(hObject, eventdata, handles);

% hFig = hObject.Parent;
% ax = findobj(hFig.Children,'type','axes');
% lines = ax(end).Children;
% for l = lines(:)'
%   l.MarkerSize = source.Value;
% end

end