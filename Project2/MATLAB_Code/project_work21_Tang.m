function A = project_work21_Tang(P1, P2, coord, elements, basis, freq)
    N_elements = length(elements);
    C = 0.3;
    
    c0 = physconst('LightSpeed');
    lambda = c0/freq;
    k_wave = 2*pi/lambda;
    mu0 = 4*pi*1e-7;
    omega = 2*pi*freq;
    eps0 = 8.854187817e-12;

    N_interior = max(max(abs(basis)));
    A = zeros(N_interior);


    for k = 1: N_elements
        for l = 1: N_elements
            node_indices_1 = elements(:, k);
            node_indices_2 = elements(:, l);
            
            p1 = coord(:, node_indices_1(1));
            p2 = coord(:, node_indices_1(2));
            p3 = coord(:, node_indices_1(3));
            c_P = (p1+p2+p3)/3;
            edge_P = [norm(p2-p1), norm(p3-p2), norm(p1-p3)];
            L_P = sum(edge_P);
            
            q1 = coord(:, node_indices_2(1));
            q2 = coord(:, node_indices_2(2));
            q3 = coord(:, node_indices_2(3));
            c_Q = (q1+q2+q3)/3;
            edge_Q = [norm(q2-q1), norm(q3-q2), norm(q1-q3)];
            L_Q = sum(edge_Q);

            if  k == l
                alok0 = alok0_cal_singular(k_wave, p1, p2, p3, P2);
            elseif norm(c_P-c_Q) < C*(L_Q+L_P)
                alok0 = alok0_cal_near_singular(k_wave, p1, p2, p3, q1, q2, q3, P2);
            else
                alok0 = alok0_cal_regular(k_wave, p1, p2, p3, q1, q2, q3, P1);
            end

            [alok1, alok2] = alok1_alok2_cal(alok0, p1, p2, p3, q1, q2, q3);

            for i1 = 1: 3

                n1 = abs(basis(i1, k));
                if n1 == 0
                    continue;
                end

                for i2 = 1: 3
                    
                    n2 = abs(basis(i2, l));
                    if n2 == 0
                        continue;
                    end

                    s1 = sign(basis(i1, k));
                    s2 = sign(basis(i2, l));
                    
                    A(n1, n2) = A(n1, n2)+s1*s2*edge_P(i1)*edge_Q(i2)*...
                        (1/(1i*omega*eps0)*alok2(i1, i2)+1i*omega*mu0*alok1(i1, i2));
                end
            end
        end
    end
end