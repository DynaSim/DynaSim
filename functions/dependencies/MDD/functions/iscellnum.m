function OUT = iscellnum(IN)
% ISCELLNUM(S) returns 1 if IN is a cell array of numerics and 0 otherwise.
    
    if iscell(IN)
        OUT = all(cellfun(@isnumeric,IN(:)));
    else
        OUT = false;
    end

end