function [matched,error_message]=CheckHostPaths(studyinfo)
%% CheckHostPaths(studyinfo)
% Purpose: compare paths on host to those set in studyinfo when batch was
% created.
% Inputs:
%   - studyinfo: DynaSim studyinfo structure
% Outputs
%   - matched: {0 or 1} (whether the paths match or not)
% 
% Paths compared:
% - path to DynaSim functions
% - path to model files

% call CheckHostPaths() in job.m before first simulation to make sure the 
% paths on the compute matchine match those in studyinfo created when the 
% study was began. 

% if paths do not match, an informative error message is added to
% .error_log (see CreateBatch() or any jobX.m script).

matched=1;
error_message='';

% locate DynaSim toolbox
dynasim_path=fileparts(fileparts(which(mfilename))); % root is one level up from directory containing this function

% locate mechanism files
[mech_paths,mech_files]=LocateModelFiles(studyinfo.base_model);

% compare DynaSim toolbox paths
if ~isequal(dynasim_path,studyinfo.paths.dynasim_functions)
  matched=0;
  error_message=sprintf('%sPath changed to DynaSim functions (expected: %s, found: %s). ',error_message,studyinfo.paths.dynasim_functions,dynasim_path);
end
  
% compare model paths
if ~isequal(unique(mech_paths),unique(studyinfo.paths.mechanisms))
  matched=0; mech_paths, studyinfo.paths.mechanisms
  error_message=sprintf('%sPath changed to model files (expected: %s, found: %s). ',error_message,[studyinfo.paths.mechanisms{:}],[mech_paths{:}]);
end  
