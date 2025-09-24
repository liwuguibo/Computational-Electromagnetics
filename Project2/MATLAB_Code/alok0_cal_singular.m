function alok0 = alok0_cal_singular(k, p1, p2, p3, N)
    [xi, eta, nu] = gausstriangle(N);

    P = [p1, p2, p3];
    Ak = 0.5*norm(cross(p2-p1, p3-p1));

    Np = length(nu);
    Nhat = [1-xi-eta; xi; eta]; 

    alok0 = zeros(3, 3);

    tol = 1e-10;
    for p = 1 : Np
        r_p = P * Nhat(:, p);
        int_vals = intRNs2(r_p, p1, p2, p3, k);

        for i = 1: 3
            N_i_p = Nhat(i, p);
            for j = 1: 3
                sum_G = 0;
                for q = 1: Np
                    r_q = P*Nhat(:, q);
                    R  = norm(r_p-r_q);

                    if R < tol
                        Gcorr = 1i * k / (4 * pi);
                    else
                        G_full = exp(1i * k * R) / (4 * pi * R);
                        Gcorr = G_full-1/(4*pi*R)+(k^2*R)/(8*pi);
                    end

                    sum_G = sum_G+nu(q)*(Nhat(j, q)*Gcorr);
                end

                contrib_p = N_i_p * (int_vals(j, 1)+int_vals(j, 2)+2*Ak*sum_G);
                alok0(i, j) = alok0(i, j)+nu(p)*contrib_p;
            end
        end
    end
    alok0 = 2*Ak*alok0;
end
