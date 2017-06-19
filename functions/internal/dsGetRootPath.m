function ds_root_path = dsGetRootPath()
  % getRootPath - get path to main ds directory
  ds_root_path = fileparts(fileparts(which('dsSimulate')));
end
