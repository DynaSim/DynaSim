function devBool = isDsDevMode()
% check if _dev directory exists in main dynasim directory

dsPath = fullfile(thisMfileDir(), '..');
devBool = isdir(fullfile(dsPath, '_dev'));

end