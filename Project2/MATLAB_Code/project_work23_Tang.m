function F = project_work23_Tang(P, XYZ, coord, elements, basis, freq, x)
    c0 = physconst('LightSpeed');
    mu0 = 4*pi*1e-7;
    lambda = c0/freq;
    k_wave = 2*pi/lambda;
    omega = 2*pi*freq;
    
    N_elements = length(elements);
    M = length(XYZ);
    
    r0 = zeros(1, M);
    
    for i = 1: M
        r0(i) = norm(XYZ(:, i));
    end
    
    r_hat = XYZ./r0;
    common_factor = -1i*omega*mu0*exp(1i*k_wave*r0)./(4*pi*r0);
    
    F_raw = zeros(1, M);
    E = zeros(3, M);
    
    for m = 1: M
        for k = 1: N_elements
    
            p1 = coord(:, elements(1, k));
            p2 = coord(:, elements(2, k));
            p3 = coord(:, elements(3, k));
            
            dlok = dlok_cal(P, p1, p2, p3, r_hat(:, m), k_wave);
            
            for j = 1: 3
    
                n = abs(basis(j, k));
                if n == 0
                    continue;
                end
    
                s_jl = sign(basis(j, k));
                idxA = elements(j, k);
                idxB = elements(mod(j, 3)+1, k);
                L_jl = norm(coord(:, idxA)-coord(:, idxB));
    
                dlok_j = dlok(:, j);
                
                far = cross(r_hat(:, m), cross(r_hat(:, m), dlok_j));
    
                c = x(n);
                
                E(:, m) = E(:, m)+common_factor(m)*(c*s_jl*L_jl)*far;
            end
        end
    
        F_raw(:, m) = norm(E(:, m));
    end
    
    F = F_raw/max(F_raw);
end