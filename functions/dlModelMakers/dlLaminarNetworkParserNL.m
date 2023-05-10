function y = dlLaminarNetworkParserNL(x, totalSize)

    c = sum(x)*4;
    z = zeros(4, 6);
    z(2:4, :) = x;
    z(1, :) = c;

    z = totalSize * z / sum(sum(z));

    y = zeros(4, 3);
    y(:, 1) = sum(z(:, 1:3), 2);
    y(:, 2) = z(:, 4);
    y(:, 3) = sum(z(:, 5:6), 2);

    y(1, :) = round(y(1, :));
    y(2:4, :) = round(y(2:4, :));
    residueN = totalSize - sum(sum(y));
    y(y < 1) = 1;

    while residueN > 0

        y(1, mod(residueN, 3)+1) = y(1, mod(residueN, 3)+1) + 1;
        residueN = residueN - 1;

    end

    while residueN < 0

        y(1, mod(residueN, 3)+1) = y(1, mod(residueN, 3)+1) - 1;
        residueN = residueN + 1;

    end

end