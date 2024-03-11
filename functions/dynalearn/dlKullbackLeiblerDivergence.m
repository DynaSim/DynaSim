function d = dlKullbackLeiblerDivergence(p, q)

    n = length(p);
    m = length(q);
    N = linspace(0, 1, n);
    M = linspace(0, 1, m);

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
    d(isnan(d) | isinf(d)) = 0;
    d = mean(d);

end