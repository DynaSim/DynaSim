function y = dlMpool1d(x, w)
    
    n = size(x, 2);
    y = x;
    for i = 1:n-w
        
        y(i) = mean(x(i:i+w));

    end
    
end