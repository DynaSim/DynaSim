function output_txt = imdpDataCursorCallback(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

%get studyID
UserData = event_obj.Target.Parent.Parent.UserData;
mdData = UserData.mdData;
keyboard

if isfield(UserData,'xLabel')
  xLabel = UserData.xLabel;
else
  xLabel = 'X';
end

if isfield(UserData,'yLabel')
  yLabel = UserData.yLabel;
else
  yLabel = 'Y';
end

if isfield(UserData,'zLabel')
  zLabel = UserData.zLabel;
else
  zLabel = 'Z';
end

pos = get(event_obj,'Position');
output_txt = {[xLabel ': ',num2str(pos(1),4)],...
    [yLabel ': ',num2str(pos(2),4)]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = [zLabel ': ',num2str(pos(3),4)];
end

xData = UserData.xData;
yData = UserData.yData;
simIdData = UserData.simIdData;

x=pos(1);
y=pos(2);
simID = simIdData( logical((xData==x) .* (yData==y)) );

output_txt{end+1} = ['simID: ', num2str(simID)];

if UserData.viewBool
  %open waveform for data point
  plotType = UserData.plotType;
  plotDir = UserData.plotDir;
  dirList = lscell(plotDir, true);
  simFiles = regexp(dirList, [plotType '.*sim' num2str(simID)]);
  pathCell = dirList(~cellfun(@isempty, simFiles));
  viewH = UserData.viewH;
  viewAxH = findobj(viewH.Children,'type','axes');
  if ~isempty(pathCell)
    filePath = fullfile(plotDir, pathCell{1});
    if exist(filePath, 'file')
%       earSuptitle(plotType, viewH)
      imshow(filePath, 'Parent', viewAxH);
    end
  else
    cla(viewAxH);
    xlim(viewAxH, [0,1])
    ylim(viewAxH, [0,1])
    text(viewAxH, 0.1,0.5,sprintf('No %s plot found for simID %i',plotType, simID))
  end
end