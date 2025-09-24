function S = project_work22_Tang(coord, elements, basis, ports, freq)
    P1 = 3;
    P2 = 4;

    A = project_work21_Tang(P1, P2, coord, elements, basis, freq);

    N_interior = max(max(abs(basis)));
    N_ports = length(ports);
    
    width = 3e-3;
    
    edge_length = zeros(N_interior, 1);
    for k = 1: length(elements)
        for i_edge = 1: 3
            p = abs(basis(i_edge, k));
            if p == 0
                continue;
            end
            
            idxA = elements(i_edge, k);
            idxB = elements(mod(i_edge, 3)+1, k);
            
            edge_length(p) = norm(coord(:, idxA) - coord(:, idxB));
        end
    end
    
    B = zeros(N_interior, N_ports);
    for m = 1: N_ports
        B(ports(m), m) = 1;
    end
    
    X = A\B;
    
    I_phys = zeros(N_interior, N_ports);
    for i = 1: N_ports
        for n = 1: N_interior
            I_phys(n, i) = X(n, i)*(edge_length(n)/2);
        end
    end
    
    Y_in = zeros(N_ports, N_ports);
    for i = 1: N_ports
      for j = 1: N_ports
        p = ports(j);
        Y_in(j,i) = -I_phys(p, i)*width;
      end
    end
    
    Z0 = 50;
    I3 = eye(N_ports);
    S = (I3-Z0*Y_in)/(I3+Z0*Y_in);
end