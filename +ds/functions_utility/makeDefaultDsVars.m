function makeDefaultDsVars()
% makeDefaultDsVars - write default DS vars to disk as text file in main DS dir.

%% Get Vars
% Get home path
if ispc
    host_name= getenv('COMPUTERNAME');
    home_path= getenv('HOMEPATH');
else
    [~, host_name]=system('echo $HOSTNAME');
    home_path= getenv('HOME');
end
host_name = strip(host_name); % remove whitespace from host_name

docs_path = fullfile(home_path, 'Documents');

% System-specific customization
dynaSimData_path = fullfile(docs_path,'DynaSimData');
demos_path = fullfile(dynaSimData_path, 'demos');

ds_path = fileparts(fileparts(which('SimulateModel')));
demos_zips_path = fullfile(ds_path, 'demos','demo_zips');

%% Write vars to disk
vars = who; % get all vars

fid = fopen(fullfile(ds_path, 'dsVars.txt'), 'w');

for thisVar = vars(:)'
  thisVar = thisVar{1};
  
  fprintf(fid, '%s = %s\r\n', thisVar,eval(thisVar));
end

fclose(fid);

end