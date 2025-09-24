function [alok3, alok4] = exercise43_Tang(p1, p2, p3, Q)

    x_p = [p1(1) p2(1) p3(1)];
    y_p = [p1(2) p2(2) p3(2)];
    
    f1 = @(x, y) 1-x-y;
    f2 = @(x, y) x;
    f3 = @(x, y) y;
    gradient_f = [-1 1 0; -1 0 1];
    
    J_matrix = [x_p(2)-x_p(1), x_p(3)-x_p(1);...
        y_p(2)-y_p(1), y_p(3)-y_p(1)];
    J = abs(det(J_matrix));
    JT_inv = inv(transpose(J_matrix));
    
    F1 = @(x, y) f1(x, y)*JT_inv*gradient_f(:, 2) - f2(x, y)*JT_inv*gradient_f(:, 1);
    F2 = @(x, y) f2(x, y)*JT_inv*gradient_f(:, 3) - f3(x, y)*JT_inv*gradient_f(:, 2);
    F3 = @(x, y) f3(x, y)*JT_inv*gradient_f(:, 1) - f1(x, y)*JT_inv*gradient_f(:, 3);
    
    [xi, eta, nu] = gausstriangle(Q);
    
    alok3 = zeros(3, 3);
    
    for i = 1: 3
        for j = 1: 3
            I = 0;
            if i == 1
                F1_local = @(x, y) F1(x, y);
            elseif i == 2
                F1_local = @(x, y) F2(x, y);
            else
                F1_local = @(x, y) F3(x, y);
            end
    
            if j == 1
               F2_local = @(x, y) F1(x, y);
            elseif j == 2
               F2_local = @(x, y) F2(x, y);
            else
               F2_local = @(x, y) F3(x, y);
            end
    
            for n = 1: length(xi)
                I = I + nu(n)*dot(F1_local(xi(n), eta(n)), F2_local(xi(n), eta(n)));
            end
            alok3(i, j) = I*J;
        end
    end
    
    alok4 = ones(3, 3)/J*2;
end