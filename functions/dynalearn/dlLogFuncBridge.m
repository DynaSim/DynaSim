function dlLogFuncBridge(logfilename)

    fileID = fopen('dlLogTempFunc.m', 'w');
    fprintf(fileID, 'function out = dlLogTempFunc(dlObj, dlArgs)\n\n\tout = %s(dlObj, dlArgs);\n\nend', logfilename);
    fclose(fileID);

end