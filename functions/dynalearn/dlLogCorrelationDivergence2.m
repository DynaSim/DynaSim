function d = dlLogCorrelationDivergence2(p, q)

    n = size(p, 1);
    m = size(p, 2);

    if isnan(p)

        d = 1e+2;
        return;

    else

        q = imresize(q, size(p));
        dlP = p / max(max(p));
        dlQ = q / max(max(q));

        % th = max(max(min(dlP, [], "all"), min(dlQ, [], "all")), 0.1);
        th = 0.1;
        dlP(dlP < th) = th;
        dlQ(dlQ < th) = th;

    end

    % figure();
    % 
    % subplot(1, 2, 1);
    % imagesc(dlP);
    % 
    % subplot(1, 2, 2);
    % imagesc(dlQ);

    d = (dlQ .* log(dlP ./ dlQ));
    d(isnan(d)) = 0;
    d(isinf(d)) = max(n, m);
    d = sum(d.^2, "all");

end