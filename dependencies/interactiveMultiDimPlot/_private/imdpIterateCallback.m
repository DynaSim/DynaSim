function imdpIterateCallback(hObject, eventdata, handles)

if hObject.Value && (~isValidFigHandle('handles.PlotPanel.figHandle') || ~handles.PlotPanel.nViewDims)
  wprintf('Cannot iterate without a Plot Panel and at least 1 ViewDim.')
  hObject.Value = 0;
  return
elseif hObject.Value% turn on
  hObject.String = 'Iterate [on]';
else % turn off
  hObject.String = 'Iterate [off]';
end

% Vars
% dimVals = handles.mdData.dimVals;
nDimVals = handles.mdData.nDimVals;

viewDims = handles.PlotPanel.viewDims;
if sum(viewDims) > 2
  viewDims(:) = 0;
end
incrDims = ~viewDims;

axDimsReversed = 1:length(incrDims);
axDimsReversed = axDimsReversed(incrDims);
axDimsReversed = flip(axDimsReversed);

% Loop
iterBool = hObject.Value;
while iterBool
  tic
  
  incrementSliders();
  
  % Check for off
  iterBool = hObject.Value;
  
  % Check for delay time
  delayTime = handles.delayBox.Value;
  iterTime = toc;
  delayTimeFinal = max( (delayTime - iterTime), 0);
  pause(delayTimeFinal)
end

  function incrementSliders()
    handles = getappdata(handles.output, 'UsedByGUIData_m');
     
    axInd = handles.PlotPanel.axInd;
    
    % check if all values arent at final index
    if ~all(axInd(incrDims) == nDimVals(incrDims))
%       TODO: recursiveIncre() instead
      for axDim = axDimsReversed
        if axInd(axDim) < nDimVals(axDim) %Then iterate
          sliderObject = handles.(handles.MainPanel.HandlesNames.sH{axDim});
          sliderObject.Value = sliderObject.Value + sliderObject.SliderStep(1)*(sliderObject.Max-sliderObject.Min);
          imdpSliderChangeCallback(sliderObject, eventdata, handles);
          break % stop incrementing others
        end
      end
    else % Reset to starting index
      % loop over slider handles and set value to min
      for iSlider = axDimsReversed
        sliderObject = handles.(handles.MainPanel.HandlesNames.sH{iSlider});
%         thisMinVal = dimVals{iSlider}(1);
%         sliderObject.Value = thisMinVal;
%         sliderObject.UserData.sibling.Value = thisMinVal;
%         sliderObject.UserData.sibling.String = num2str(thisMinVal);
        
        sliderObject.Value = -inf;
        handles = imdpSliderChangeCallback(sliderObject, eventdata, handles);
        
        % Update handles structure
        guidata(hObject, handles);
      end
      
      % replot
      handles = imdpPlot(hObject, eventdata, handles);
      
      % Update handles structure
      guidata(hObject, handles);
    end
  end %incrementSliders

end
