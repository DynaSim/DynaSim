function devBool = dsIsDevMode()
% check if _dev directory exists in main dynasim directory

dsGetRootPath = fullfile(thisMfileDir(), '..');
devBool = isfolder(fullfile(dsGetRootPath, '_dev'));

end
