function d = dlLogCorrelationDivergence(p, q)

    n = length(p);
    m = length(q);
    N = linspace(0, 1, n);
    M = linspace(0, 1, m);

    if isnan(p)

        d = 1e+2;
        return;

    end

    if n > m

        dlQ = interp1(M, q, N);
        dlP = p / max(max(p));
        dlQ = dlQ / max(max(dlQ));

    elseif m > n

        dlP = interp1(N, p, M);
        dlP = dlP / max(max(dlP));
        dlQ = q / max(max(q));

    else

        dlP = p / max(max(p));
        dlQ = q / max(max(q));

    end

    d = abs(log(dlP ./ dlQ));
    d(isnan(d)) = 0;
    d(isinf(d)) = max(n, m);
    d = mean(d);

end