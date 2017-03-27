

function out = getpath(query)
% out = getpath(query)
% Author David Stanley stanleyd@bu.edu
% Returns the value of the path variable specificed in query

    mypaths
    eval(['out = ' query ';']);
    
end

