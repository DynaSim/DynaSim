function ds_root_path = getRootPath()
  % getRootPath - get path to main ds directory
  ds_root_path = fileparts(fileparts(which('ds.simulateModel')));
end
