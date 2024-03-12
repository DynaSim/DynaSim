function d = dlKullbackLeiblerDivergence(p, q)

    n = length(p);
    m = length(q);
    N = linspace(0, 1, n);
    M = linspace(0, 1, m);

    if isnan(p) || isinf(p)

        d = 1e+6;
        return;

    end

    if n > m

        dlQ = interp1(M, q, N);
        dlP = p / sum(p);
        dlQ = dlQ / sum(dlQ);

    elseif m > n

        dlP = interp1(N, p, M);
        dlP = dlP / sum(dlP);
        dlQ = q / sum(q);

    else

        dlP = p / sum(p);
        dlQ = q / sum(q);

    end

    d = dlP .* log(dlP ./ dlQ);
    d(isnan(d)) = 0;
    d(isinf(d)) = max(n, m);
    d = mean(d);

end