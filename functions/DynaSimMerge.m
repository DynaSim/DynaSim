
function data_merged = DynaSimMerge(data1,data2)
    % Merge two Dynasim structures created from different simulations
    
    xp1 = All2xPlt(data1);
    xp2 = All2xPlt(data2);
    
    xp1 = xp1.merge(xp2);
    
    if isfield(data1,'plot_files') && isfield(data2,'plot_files')
        data_merged = xPlt2DynaSimImage(xp1);
    elseif ~isfield(data1,'plot_files') && ~isfield(data2,'plot_files')
        data_merged = xPlt2DynaSim(xp1);
    else
        error('Unknown input types');
    end
    
end

