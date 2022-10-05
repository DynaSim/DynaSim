function [y, l] = dlAccuracyBastos2020Task(dlObj, args)

    x = dlObj.dlLastOutputs;
    x = cell2mat(x(1:3));
    l = max(x);
    y = find(x == l);

end

