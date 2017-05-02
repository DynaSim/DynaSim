function varOutput = getConfig(query)
% GETCONFIG - returns the value of the variable specificed as string in argument
%
% Usage: pathOutput = getConfig(pathQueryString)

  dsVarsFile = fullfile(ds.getRootPath(),'dsConfig.txt');

  if ~exist(dsVarsFile, 'file')
    ds.makeDefaultConfig(); % bring all path vars into namespace
  end
  
  fid = fopen(dsVarsFile);
  dsVars = textscan(fid, '%s = %s');
  fclose(fid);
  
  varCell = dsVars{2}(~cellfun(@isempty,strfind(dsVars{1}, query)));
  
  if ~isempty(varCell)
    varOutput = varCell{1}; % eval variable string as variable
  else
    varOutput = [];
    if isempty(varOutput); warning('Requested path not found. dsVars.txt is possibly corrupt. Try deleting dsConfig.txt and running ds.makeDefaultConfig()');
  end
    
end
