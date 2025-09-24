function [fc, xc] = exercise34_Tang(coord, elements, nodes)

x_p = coord(1, :);
y_p = coord(2, :);

N_triangle = length(elements);
N_node = length(nodes);

[xi, eta, nu] = gausstriangle(3);

A = zeros(N_node, N_node);
B = zeros(N_node, N_node);

f = @(x,y) [1-x-y, x, y];
g = [-1, 1, 0; -1, 0, 1];

for k = 1: N_triangle
    alok1 = zeros(3, 3);
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
            I_alok1 = 0;
            I_alok2 = 0;
            for n = 1: length(xi)
                temp = f(xi(n), eta(n));
                I_alok1 = I_alok1 + nu(n)*temp(i)*temp(j);

                I_alok2 = I_alok2 + nu(n)*dot((JT_inv*g(:, i)), JT_inv*g(:, j));
            end
            alok1(i, j) = I_alok1*J;
            alok2(i, j) = I_alok2*J;
        end
    end

    for i = 1: 3
        for j = 1: 3
            A(triangle_vertices(i), triangle_vertices(j)) = ...
                A(triangle_vertices(i), triangle_vertices(j)) + alok2(i, j);

            B(triangle_vertices(i), triangle_vertices(j)) = ...
                B(triangle_vertices(i), triangle_vertices(j)) + alok1(i, j);
        end
    end
end

[V, D] = eig(A, B);

D_length = length(D);
d = zeros(1, D_length);
for i = 1: D_length
    d(i) = D(i, i);
end

kc = zeros(1, 6);
fc = zeros(1, 6);
xc = zeros(N_node, 6);
c0 = physconst('LightSpeed');

for i = 1: 6
    kc(i) = sqrt(D(i+1, i+1));
    fc(i) = c0*kc(i)/(2*pi);
    xc(:,i) = V(:, i+1);
end