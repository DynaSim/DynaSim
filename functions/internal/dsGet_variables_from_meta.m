function varlabels = dsGet_variables_from_meta(xp)
    varlabels = {};
    if ~isempty(xp.meta.dynasim)
        if ~isempty(xp.meta.dynasim.labels)
            myabels = xp.meta.dynasim.labels;
            mylabels = myabels(~strcmpi(myabels,'time')); % Remove time
            inds = cellfun(@(s) strfind(s,'_'),mylabels,'UniformOutput',0);
            varlabels = cellfun(@(s,ind) s(ind+1:end),mylabels,inds,'UniformOutput',0);
        end
    end
end
