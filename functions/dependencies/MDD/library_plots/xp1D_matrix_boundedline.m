

function hxp = xp1D_matrix_boundedline (xp, op)
    % xp must be 1xN (e.g. 1 dimensional)
    
    hxp = struct;
    
    if nargin < 2
        op = struct;
    end
    
    if isempty(op); op = struct; end;
    
    op = struct_addDef(op,'xlims',[]);
    op = struct_addDef(op,'ylims',[]);
    op = struct_addDef(op,'meanfunc',@(x) mean(x,2));
    op = struct_addDef(op,'errfunc',@(x) std(x,[],2) ./ (sqrt(size(x,2)) * ones(size(x,1),1)));
    
    xlims = op.xlims;
    ylims = op.ylims;
    meanfunc = op.meanfunc;
    errfunc = op.errfunc;
    
    N = length(xp.data);
    
    % Get 1st data point for estimating size
    i=1;
    if ~isempty(xp.meta.datainfo(1).values)
        t = xp.meta.datainfo(1).values;
    else    
        t = 1:length(xp.data{i});
    end
    t = double(t);
    muarr = zeros(length(t),N);
    errarr = zeros(length(t),1,N);
    
    for i = 1:N
        if ~isempty(xp.meta.datainfo(1).values)
            t = xp.meta.datainfo(1).values;
        else    
            t = 1:length(xp.data{i});
        end
        t = double(t);
        d = double(xp.data{i});
        
        
        if ~isempty(d)
            if ismatrix(d)
                mu = meanfunc(d); mu=mu(:);
                err = errfunc(d); err=err(:);
                muarr(:,i) = mu;
                errarr(:,1,i) = err;
                hold on;
                %hxp.hcurr = boundedline(t,mu,[err err],'alpha');
                hxp.hcurr = plot(t,mu,'LineWidth',2);
            else
                error('Too many dimensions');
            end
        end
    end
    
    errarr = repmat(errarr,[1,2,1]);
    hxp.hcurrErr = boundedline(repmat(t(:),[1,N]),muarr,errarr,'alpha');
    
    if ~isempty(xlims); xlim(xlims); end
    if ~isempty(ylims); ylim(ylims); end
    
end


