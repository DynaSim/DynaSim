function ds_path = dsPath()
  ds_path = fileparts(fileparts(which('SimulateModel')));
end