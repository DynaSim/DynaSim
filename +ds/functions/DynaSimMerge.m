
function data_merged = DynaSimMerge(data1,data2)
    % Merge two Dynasim structures created from different simulations
    
    xp1 = all2xPlt(data1);
    xp2 = all2xPlt(data2);
    
    xp_merged = merge(xp1,xp2);
    xp_merged.meta = xp1.meta; % Warning - need to make sure metadata merges properly; not yet implemented
    
    % Sort everything except populations and variables
    inds = true(1,ndims(xp_merged));
    inds(xp_merged.findaxis('populations')) = false;
    inds(xp_merged.findaxis('variables')) = false;
    xp_merged = xp_merged.sortAxis(find(inds));
    
    if isfield(data1,'plot_files') && isfield(data2,'plot_files')
        data_merged = xPlt2DynaSimImage(xp_merged);
    elseif ~isfield(data1,'plot_files') && ~isfield(data2,'plot_files')
        data_merged = xPlt2DynaSim(xp_merged);
    else
        error('Unknown input types');
    end
    
end

