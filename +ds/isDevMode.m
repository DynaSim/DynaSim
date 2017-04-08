function devBool = isDevMode()
% check if _dev directory exists in main dynasim directory

ds.getRootPath = fullfile(thisMfileDir(), '..');
devBool = isdir(fullfile(ds.getRootPath, '_dev'));

end