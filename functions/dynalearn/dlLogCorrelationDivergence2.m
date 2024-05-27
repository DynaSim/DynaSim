function d = dlLogCorrelationDivergence2(p, q)

    n = size(p, 1);
    m = size(p, 2);

    if isnan(p)

        d = 1e+2;
        return;

    else

        q = imresize(q, size(p));
        dlP = p / sum(mean(p, 2));
        dlQ = q / sum(mean(q, 2));

        th = max(min(dlP, [], "all"), min(dlQ, [], "all"));
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
    d = sum(d, "all")^2;

end