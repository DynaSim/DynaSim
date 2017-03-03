function varargout = imdpSliderChangeCallback(hObject, eventdata, handles)
% Fn sets both slider and edit box to nearest value in data

% Slider handle number
hNum = regexp(hObject.Tag, '(\d+)', 'tokens');
hNum = str2double(hNum{1});

%% Find Nearest Value to New Value
lastVal = hObject.UserData.lastVal;
switch hObject.Style
  case 'slider'
    if ~isnumeric(eventdata)
      newVal = hObject.Value;
    else %from scroll
      newVal = lastVal + eventdata * (hObject.SliderStep(1))*(hObject.Max-hObject.Min);
    end
  case 'edit'
    newVal = str2double(hObject.String);
end
% varInd = hObject.UserData.varInd;
% varName = hObject.UserData.varName;
thisDimVals = handles.mdData.dimVals{hNum};
nVals = length(thisDimVals);
% nearestVal = handles.data.Table.(varName)(nearest(handles.data.Table.(varName), newVal));
nearestValInd = nearest(thisDimVals, newVal);
nearestVal = thisDimVals(nearestValInd);


% increment if newVal same as lastVal
if nearestVal == lastVal
  changeSign = sign(newVal - lastVal);
%   lastIndNearestVal = find(handles.data.Table.(varName)==nearestVal, 1, 'last');
%   nearestVal = handles.data.Table.(varName)(lastIndNearestVal + changeSign);

  if ~changeSign % == 0
    return % same val so do nothing and return
  end

  % update nearest
  nearestValInd = nearestValInd + changeSign;
  
  if nearestValInd > nVals % correct overflow
    nearestValInd = nVals;
  elseif nearestValInd < 1 % correct underflow
    nearestValInd = 1;
  end

  nearestVal = thisDimVals(nearestValInd);
end

%% Change Slider Vals
switch hObject.Style
  case 'slider'
    % change slider to closest value
    hObject.Value = nearestVal;
    hObject.UserData.lastVal = nearestVal;
    
    % change slider val
    hObject.UserData.sibling.Value = nearestVal;
    hObject.UserData.sibling.String = num2str(nearestVal);
    hObject.UserData.sibling.UserData.lastVal = nearestVal;
  case 'edit'
    % change slider to closest value
    hObject.UserData.sibling.Value = nearestVal;
    
    % change slider val
    hObject.Value = nearestVal;
    hObject.String = num2str(nearestVal);
end

% Add index
handles.PlotPanel.axInd(hNum) = nearestValInd;

%% Replot
% Check if reset from iterate or hit boundaries
if newVal > -inf && nearestVal ~= lastVal
  imdpPlot(hObject, eventdata, handles);
end

% Update handles structure
guidata(hObject, handles);

if nargout > 0
  varargout{1} = handles;
end

end