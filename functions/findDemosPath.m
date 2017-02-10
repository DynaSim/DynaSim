function path2demosDir = findDemosPath
  path2demosFile = which('demos');
  path2demosDir = fileparts(path2demosFile);
end