function alok0 = exercise84_Tang(k, p1, p2, p3, q1, q2, q3, xi, eta, nu)

    P = [p1, p2, p3];
    Q = [q1, q2, q3];

    Ak = 0.5*norm(cross(p2-p1, p3-p1));
    Al = 0.5*norm(cross(q2-q1, q3-q1));
    
    Np = length(nu);
    Nhat = [1 - xi - eta; xi; eta];

    alok0 = zeros(3, 3);

    for p = 1: Np
        rk = P*Nhat(:, p);

        intvals = intRNs2(rk, q1, q2, q3, k);

        for i = 1: 3
            for j = 1: 3
                sum_G = 0;

                for q = 1: Np
                    rl = Q*Nhat(:, q);
                    R = norm(rk-rl);

                    G = exp(1i*k*R)/(4*pi*R);
                    G_corr = G - 1/(4*pi*R) + (k^2*R)/(8*pi);

                    sum_G = sum_G + nu(q)*G_corr*Nhat(j, q);
                end
                alok0(i, j) = alok0(i, j) + ...
                    nu(p)*Nhat(i, p)*(intvals(j, 1)+intvals(j, 2)+2*Al*sum_G);
            end
        end
    end
    
    alok0 = 2*Ak*alok0;
end
