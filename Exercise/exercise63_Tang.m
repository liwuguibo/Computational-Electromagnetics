function x =  exercise63_Tang(coord, elements)
    N_triangle = length(elements);
    A = zeros(N_triangle, N_triangle);
    b = zeros(N_triangle, 1);

    NTa = N_triangle / 2;

    eps0 = 8.854187817e-12;
    k = 0;

    % Precompute centers and areas
    centers = zeros(3, N_triangle);
    areas = zeros(1, N_triangle);

    for m = 1: N_triangle
        verts = coord(:, elements(:, m));
        centers(:, m) = mean(verts, 2);
        v1 = verts(:, 2) - verts(:, 1);
        v2 = verts(:, 3) - verts(:, 1);
        areas(m) = norm(cross(v1, v2))/2;
    end
    
    for m = 1:N_triangle
        for n = 1:N_triangle
            if m == n
                p1 = coord(:, elements(1,m));
                p2 = coord(:, elements(2,m));
                p3 = coord(:, elements(3,m));
                int_val = intRNs2(centers(:,m), p1, p2, p3, k);
                A(m,m) = (1/eps0) * areas(m) * sum(int_val(:,1));
            else
                G = 1 / (4 * pi * norm(centers(:,m) - centers(:,n)));
                A(m,n) = (1/eps0) * areas(m) * areas(n) * G;
            end
        end
    end

    % Assemble b vector
    for m = 1: N_triangle
        if m <= NTa
            b(m) = 5 * areas(m);    % Inner conductor potential
        else
            b(m) = 0;               % Outer conductor potential
        end
    end

    x = A\b;
end