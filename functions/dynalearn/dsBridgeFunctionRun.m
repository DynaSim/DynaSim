function dsBridgeFunctionRun(dlOutputType)

    tempfuncname = 'dlTempFuncBridge.m';
    fileID = fopen(tempfuncname, 'w');
    fprintf(fileID, 'function out = dlTempFuncBridge(dlOutputs)\n\n\tout = %s(dlOutputs);\n\nend', dlOutputType);
    fclose(fileID);

end