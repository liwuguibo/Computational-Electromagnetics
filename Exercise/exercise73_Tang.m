function F = exercise73_Tang(freq, L, P ,XYZ)
    % Param Definition
    c0 = physconst('LightSpeed');
    lambda = c0/freq;
    k = 2*pi/lambda;

    [xi, nu] = gausslegendre(P);

    z_nodes = (L/2)*xi;
    
    % Field points
    M = length(XYZ);
    F_raw = zeros(1, M);

    for m = 1: M
        r_vec = XYZ(:, m);
        r = norm(r_vec);
        
        r_hat = r_vec/r;
        
        phase = exp(-1i*k*z_nodes*r_hat(3));
        
        inner_sum = sum(nu .* phase .* sin(pi*xi));
        
        z_hat = [0; 0; 1];
        
        E_dir = cross(r_hat, cross(r_hat, z_hat));
        E_vec = E_dir * L * inner_sum;
        F_raw(m) = norm(E_vec);
    end
    
    F = F_raw/max(F_raw);
end