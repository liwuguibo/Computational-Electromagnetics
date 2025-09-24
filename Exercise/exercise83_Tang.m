function alok0 = exercise83_Tang(k, p1, p2, p3, q1, q2, q3, xi, eta, nu)

    P = [p1, p2, p3];
    Q = [q1, q2, q3];

    Ap = 0.5*norm(cross(p2-p1, p3-p1));
    Aq = 0.5*norm(cross(q2-q1, q3-q1));

    Np = length(nu); 
    
    Nhat = [1 - xi - eta; xi; eta];

    alok0 = zeros(3, 3);

    for p = 1: Np

        rk = P*Nhat(:, p);

        for q = 1: Np

            rl = Q*Nhat(:, q);
            
            R = norm(rk-rl);
            G = exp(1i*k*R)/(4*pi*R);
            
            for i = 1:3
                for j = 1:3
                    alok0(i,j) = alok0(i,j) + ...
                        nu(p)*nu(q)*Nhat(i, p)*G*Nhat(j, q);
                end
            end
        end
    end
    alok0 = 4*Ap*Aq*alok0;
end