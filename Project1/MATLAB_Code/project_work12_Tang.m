function [beta, fc, xET, xEZ] = project_work12_Tang(A, B, C, D, edges, nodes, f0)

c0 = physconst('LightSpeed');
lambda0 = c0/f0;
k0 = 2*pi/lambda0;

int_edges = find(edges == 0);
int_nodes = find(nodes == 0);

Aii = A(int_edges, int_edges);

Bii = sparse(B(int_edges, int_edges));
Cii = sparse(C(int_edges, int_nodes));
Dii = sparse(D(int_nodes, int_nodes));

M = Bii - Cii * (Dii \ Cii.');

M = full(M);

[eigV, eigD] = eig(Aii, M);

gamma_all = diag(eigD);

[gamma_sorted, idx] = sort(gamma_all, 'ascend', 'ComparisonMethod','real');
gamma_selected = gamma_sorted(1:6);
beta = sqrt(-gamma_selected);
fc = sqrt(k0^2 - beta.^2)*c0/(2*pi);

xET_sel = eigV(:, idx(1:6));

xEZ_sel = zeros(length(int_nodes), 6);

for i = 1:6
    xEZ_sel(:, i) = 1i * beta(i) * (Dii \ (Cii.' * xET_sel(:, i)));
end

xET = zeros(length(edges), 6);
xET(int_edges, :) = xET_sel;
for i = 1: 6
    xET_max = max(abs(xET(:, i)));
    xET(:, i) = xET(:, i)/xET_max;
end

xEZ = zeros(length(nodes), 6);
xEZ(int_nodes, :) = xEZ_sel;
for i = 1: 6
    xEZ_max = max(abs(xEZ(:, i)));
    xEZ(:, i) = xEZ(:, i)/xEZ_max;
end