function dlLogFuncBridge(logfilename)

    fileID = fopen('dlLogTempFunc.m', 'w');
    fprintf(fileID, 'function [out1, out2] = dlLogTempFunc(dlObj, dlArgs)\n\n\t[out1, out2] = %s(dlObj, dlArgs);\n\nend', logfilename);
    fclose(fileID);

end