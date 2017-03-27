function varOutput = getDsVar(query)
% GETDSVAR - returns the value of the variable specificed as string in argument
%
% Usage: pathOutput = getPath(pathQueryString)

  dsVarsFile = fullfile(dsPath(),'dsVars.txt');

  if ~exist(dsVarsFile, 'dir')
    makeDefaultDsVars(); % bring all path vars into namespace
  end
  
  fid = fopen(dsVarsFile);
  dsVars = textscan(fid, '%s = %s');
  fclose(fid);
  
  varCell = dsVars{2}(~cellfun(@isempty,strfind(dsVars{1}, query)));
  
  if ~isempty(varCell)
    varOutput = varCell{1}; % eval variable string as variable
  else
    varOutput = [];
  end
    
end
