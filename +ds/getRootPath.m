function ds_path = getRootPath()
  % getRootPath - get path to main ds directory
  ds_path = fileparts(fileparts(which('ds.simulateModel')));
end
