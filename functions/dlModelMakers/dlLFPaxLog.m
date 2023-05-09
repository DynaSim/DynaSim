function [y, l] = dlLFPaxLog(dlObj, args)

    dlPotentialIndices = contains(dlObj.dlVariables, '_V');
    dlPotentialIndices(1) = 1;
    dlPotentials = cell2mat(dlObj.dlOutputs(dlPotentialIndices));
    dlPotentials = dlPotentials(2:10:end, :);
    
    y = [dlObj.dlOutputLog, dlPotentials];
    l = [];
    
end
