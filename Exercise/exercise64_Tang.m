function phi = exercise64_Tang(xyz, coord, elements, x)
    N_triangle = length(elements);
    N_fields = length(xyz);

    phi = zeros(N_fields, 1);

    eps0 = 8.854187817e-12;
    coef = 1 / (4 * pi * eps0);

    centers = zeros(3, N_triangle);
    areas = zeros(1, N_triangle);
    
    for n = 1: N_triangle
        verts = coord(:, elements(:, n));
        centers(:, n) = mean(verts, 2);
        v1 = verts(:,2) - verts(:,1);
        v2 = verts(:,3) - verts(:,1);
        areas(n) = norm(cross(v1, v2)) / 2;
    end

    for p = 1: N_fields
        r = xyz(:, p);
        acc = 0;
        for n = 1: N_triangle
            R = norm(r - centers(:, n));
            acc = acc + x(n) * areas(n) / R;
        end
        phi(p) =  acc*coef;
    end
end