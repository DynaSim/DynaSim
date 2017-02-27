function path2demosDir = findDemosPath
%FINDDEMOSPATH - Finds first instance of demos location on MATLAB path
%
% TODO: only works if DynaSim is highest on the path...

  path2demosFile = which('demos');
  path2demosDir = fileparts(path2demosFile);
end
