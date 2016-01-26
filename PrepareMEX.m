function mexfile=PrepareMEX(file)
% Purpose: take a solver m-file and compile it using the Matlab coder.
% See also: GetSolveFile

% Create a MEX configuration object
cfg = coder.config('mex');
% Turn on dynamic memory allocation
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
% Generate MEX function
eval(['codegen -d codemex -config cfg ',sprintf(file)]);
