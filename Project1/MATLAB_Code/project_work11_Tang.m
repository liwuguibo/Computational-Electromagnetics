function [A, B, C, D] = project_work11_Tang(coord, elements, basis, edges, nodes, f0)

c0 = physconst('LightSpeed');
lambda0 = c0/f0;
k0 = 2*pi/lambda0;

N_edges = length(edges);
N_nodes = length(coord);
N_triangle = length(elements);

Q = 3;

A = zeros(N_edges, N_edges);
B = zeros(N_edges, N_edges);
C = zeros(N_edges, N_nodes);
D = zeros(N_nodes, N_nodes);

for k = 1: N_triangle
    
    m = elements(:, k);
    b = basis(:, k);
    s = sign(b);
    n = abs(b);
    
    p1 = coord(:, m(1));
    p2 = coord(:, m(2));
    p3 = coord(:, m(3));

    [alok1, alok2, alok3, alok4, alok6] = calc_aloks_Tang(p1, p2, p3 ,Q);


    
    for i = 1: 3
        for j = 1: 3

            if edges(n(i)) == 0 && edges(n(j)) == 0
                A(n(i), n(j)) = A(n(i), n(j)) + s(i)*s(j)*alok4(i, j) - s(i)*s(j)*k0^2*alok3(i, j);
                B(n(i), n(j)) = B(n(i), n(j)) + s(i)*s(j)*alok3(i, j);
            end

            if edges(n(i)) == 0 && nodes(m(j)) == 0
                C(n(i), m(j)) = C(n(i), m(j)) + s(i)*alok6(i, j);
            end

            if nodes(m(i)) == 0 && nodes(m(j)) == 0
                D(m(i), m(j)) = D(m(i), m(j)) + alok2(i, j) - k0^2*alok1(i, j);
            end
        end
    end
end