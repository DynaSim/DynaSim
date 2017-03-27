

function out = getpath(query)
% out = getpath(query)
% Returns the value of the path variable specificed in query

    mypaths
    eval(['out = ' query ';']);
    
end

