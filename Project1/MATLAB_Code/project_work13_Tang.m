function [E, H, xy] = project_work13_Tang(xETj, xEzj, coord, elements, basis, f0, betaj)

omega = 2*pi*f0;
mu0 = 4 * pi * 1e-7;

N_triangle = length(elements); 
E = zeros(3, N_triangle);
H = zeros(3, N_triangle);
xy = zeros(2, N_triangle);

f1 = @(x, y) 1-x-y;
f2 = @(x, y) x;
f3 = @(x, y) y;
gradient_f = [-1 1 0; -1 0 1];

xi = 1/3;
eta = 1/3;

for k = 1: N_triangle
    
    coord_indice = elements(:, k);
    p = coord(:, coord_indice);
    edges_local = basis(:,k);
    s = sign(edges_local);
    eids = abs(edges_local);

    xy(:, k) = mean(p, 2);
    x_p = p(1, :);
    y_p = p(2, :);

    J_matrix = [x_p(2)-x_p(1), x_p(3)-x_p(1); y_p(2)-y_p(1), y_p(3)-y_p(1)];
    J = abs(det(J_matrix));
    JT_inv = inv(transpose(J_matrix));

    F1 = @(x, y) f1(x, y)*JT_inv*gradient_f(:, 2) - f2(x, y)*JT_inv*gradient_f(:, 1);
    F2 = @(x, y) f2(x, y)*JT_inv*gradient_f(:, 3) - f3(x, y)*JT_inv*gradient_f(:, 2);
    F3 = @(x, y) f3(x, y)*JT_inv*gradient_f(:, 1) - f1(x, y)*JT_inv*gradient_f(:, 3);
    
    ET = s(1)*xETj(eids(1)) * F1(xi, eta) + ...
         s(2)*xETj(eids(2)) * F2(xi, eta) + ...
         s(3)*xETj(eids(3)) * F3(xi, eta);
    
    EZ = sum(xEzj(coord_indice))/3;

    E(:, k) = [ET; EZ];
    
    gradN1 = JT_inv * gradient_f(:,1);
    gradN2 = JT_inv * gradient_f(:,2);
    gradN3 = JT_inv * gradient_f(:,3);
    
    HZ = (gradN2(1)*ET(2) - gradN1(2)*ET(1))/(1j*omega*mu0);

    ez_cross_ET = [-ET(2); ET(1)];
    
    gradEz = xEzj(coord_indice(1)) * gradN1 + ...
             xEzj(coord_indice(2)) * gradN2 + ...
             xEzj(coord_indice(3)) * gradN3;

    ez_cross_gradEz = [-gradEz(2); gradEz(1)];

    HT = (1i*betaj*ez_cross_ET - ez_cross_gradEz)/(1i*omega*mu0);

    H(:, k) = [HT; HZ]; 
end