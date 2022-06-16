function [B1, B2, B3, T1, T2, T3, TrB, TrT] = dlTrialBlockGenerator(dlInputParameters, dlTargetParameters, n, m)

    % Only based on DynaLearn*
    % dlInputParameters: all different input parameters
    % dlTargetParameters : all target parameters based on inputs
    % n = 50, m = 50 size of block and trial sessions
    
    B1 = cell(1, n);B2 = cell(1, n);B3 = cell(1, n);
    T1 = cell(1, n);T2 = cell(1, n);T3 = cell(1, n);
    TrB = cell(1, m);TrT = cell(1, m);
    
    for i = 1:n
        
        B1{i} = dlInputParameters{1};
        B2{i} = dlInputParameters{2};
        B3{i} = dlInputParameters{3};
        
        T1{i} = dlTargetParameters{1};
        T2{i} = dlTargetParameters{2};
        T3{i} = dlTargetParameters{3};
        
    end
    
    for i = 1:m
        
        k = randi(3);
        TrB{i} = dlInputParameters{k};
        TrT{i} = dlTargetParameters{k};
        
    end
    
end