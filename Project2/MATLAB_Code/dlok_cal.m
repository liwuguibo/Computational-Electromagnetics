function dlok = dlok_cal(P, p1, p2, p3, r, k_wave)
    r_hat = r/norm(r);

    v12 = p2-p1;
    v13 = p3-p1;
    n_vec = cross(v12, v13);
    area2 = norm(n_vec);
    n_hat = n_vec/norm(n_vec);

    gradN = zeros(3,3);
    gradN(:, 1) = cross(p2-p3, n_hat)/area2;
    gradN(:, 2) = cross(p3-p1, n_hat)/area2;
    gradN(:, 3) = cross(p1-p2, n_hat)/area2;
    
    [xi, eta, nu] = gausstriangle(P);

    Np = length(nu);

    Rprim = zeros(3, Np);

    for pp = 1: Np
        Rprim(:, pp) = p1*(1-xi(pp)-eta(pp))+p2*xi(pp)+p3*eta(pp);
    end

    dlok = zeros(3, 3);

    for j_loc = 1: 3
        
        iN = j_loc;
        kN = mod(j_loc, 3)+1;
        
        dlok_j = zeros(3,1);

        for pp = 1: Np
            rprime = Rprim(:, pp);

            phase = exp(-1i*k_wave*(r_hat.*rprime));

            Npt = [1 - xi(pp) - eta(pp); xi(pp); eta(pp)];

            vjl = -cross(n_hat, (Npt(iN).*gradN(:,kN)-Npt(kN).*gradN(:,iN)));

            dlok_j = dlok_j+(area2)*nu(pp)*phase.*vjl;
        end
        dlok(:, j_loc) = dlok_j;
    end
end