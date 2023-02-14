function y = dlTrialBlockGenerator(dlInputParameters, dlTargetParameters, n, m)

    % Only based on DynaLearn*
    % dlInputParameters: all different input parameters
    % dlTargetParameters : all target parameters based on inputs
    % n = 50, m = 50 size of block and trial sessions
    
    y = struct();
    y.B1 = cell(1, n);y.B2 = cell(1, n);y.B3 = cell(1, n);
    y.T1 = cell(1, n);y.T2 = cell(1, n);y.T3 = cell(1, n);
    y.TrB = cell(1, m);y.TrT = cell(1, m);
    
    for i = 1:n
        
        y.B1{i} = dlInputParameters{1};
        y.B2{i} = dlInputParameters{2};
        y.B3{i} = dlInputParameters{3};
        
        y.T1{i} = dlTargetParameters{1};
        y.T2{i} = dlTargetParameters{2};
        y.T3{i} = dlTargetParameters{3};
        
    end
    
    for i = 1:m
        
        k = randi(3);
        y.TrB{i} = dlInputParameters{k};
        y.TrT{i} = dlTargetParameters{4};
        
    end
    
end