
function data_merged = mergeData(data1,data2)
    % MERGEDATA - Merge two Dynasim structures created from different simulations
    
    xp1 = ds.all2MDD(data1);
    xp2 = ds.all2MDD(data2);
    
    xp_merged = mergeData(xp1,xp2);
    xp_merged.meta = xp1.meta; % Warning - need to make sure metadata merges properly; not yet implemented
    
    % Sort everything except populations and variables
    inds = true(1,ndims(xp_merged));
    inds(xp_merged.findaxis('populations')) = false;
    inds(xp_merged.findaxis('variables')) = false;
    xp_merged = xp_merged.sortAxis(find(inds));
    
    if isfield(data1,'plot_files') && isfield(data2,'plot_files')
        data_merged = MDD2dsImage(xp_merged);
    elseif ~isfield(data1,'plot_files') && ~isfield(data2,'plot_files')
        data_merged = MDD2ds(xp_merged);
    else
        error('Unknown input types');
    end
    
end
