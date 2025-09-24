function x = exercise33_Tang(coord, elements, nodes, U, epsr)

x_p = coord(1, :);
y_p = coord(2, :);

epsilon_r = epsr;
U_outer = U;
U_inner = 0;

N_triangle = length(elements);
N_node = length(nodes);

[xi, ~, nu] = gausstriangle(3);

g = [-1, 1, 0; -1, 0, 1];

A = zeros(N_node, N_node);
b = zeros(N_node, 1);

for k = 1: N_triangle
    alok2 = zeros(3, 3);

    triangle_vertices = elements(:, k);

    J_matrix = [x_p(triangle_vertices(2)) - x_p(triangle_vertices(1)),...
        x_p(triangle_vertices(3)) - x_p(triangle_vertices(1));...
        y_p(triangle_vertices(2)) - y_p(triangle_vertices(1)),...
        y_p(triangle_vertices(3)) - y_p(triangle_vertices(1))];
    J = abs(det(J_matrix));
    JT_inv = inv(transpose(J_matrix));
    
    for i = 1: 3
        for j = 1: 3
            I_alok2 = 0;
            for n = 1: length(xi)
                I_alok2 = I_alok2 + nu(n)*dot((JT_inv*g(:, i)), JT_inv*g(:, j));
            end
            alok2(i, j) = I_alok2*J;
        end
    end
    for i = 1: 3
        for j = 1: 3
            A(triangle_vertices(i), triangle_vertices(j)) = ...
                A(triangle_vertices(i), triangle_vertices(j)) + epsilon_r*alok2(i, j);
        end
    end
end

for m = 1: N_node
    if nodes(m) == 1
        for n = 1: N_node
            A(m, n) = 0;
        end
        A(m, m) = 1;
        b(m) = U_outer;
    elseif nodes(m) == 2
        for n = 1: N_node
            A(m, n) = 0;
        end
        A(m, m) = 1;
        b(m) = U_inner;
    end
end
x = A\b;