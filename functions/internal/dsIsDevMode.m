function devBool = dsIsDevMode()
% check if _dev directory exists in main dynasim directory

dsGetRootPath = fullfile(thisMfileDir(), '..');
devBool = isdir(fullfile(dsGetRootPath, '_dev'));

end