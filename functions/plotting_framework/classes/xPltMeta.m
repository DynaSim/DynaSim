
classdef xPltMeta
    % Metadata information. This is inserted into the "meta" property of
    % xPlt class

    properties
        axis_name       % (optional) 1x1 - string naming the axis in question (e.g. time, voltage, cell number, etc)
        axis_values     % (optional) 1xN - can be numerical matrix or cell array of strings describing axis in question
        axis_struct     % (optional) 1xN structure array of additional information
    end
    
    methods
        function xp_meta = inds(xp_meta,ind)
            xp_meta.axis_values = xp_meta.axis_values(ind,:);
            xp_meta.axis_text = xp_meta.axis_text(ind);
            xp_meta.axis_metadata = xp_meta.axis_metadata(ind);
        end
        
    end
end
