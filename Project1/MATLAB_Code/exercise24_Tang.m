function [alok1, alok2] = exercise24_Tang(p1, p2, p3, Q)

    % Get the position of the nodes
    x = [p1(1), p2(1), p3(1)];
    y = [p1(2), p2(2), p3(2)];
    
    % Compute integration points and weights
    [xi, eta, nu] = gausstriangle(Q);
    
    % Define integral function
    f = @(x,y) [1-x-y, x, y];
    g = [-1, 1, 0; -1, 0, 1];
    
    % Compute Jacobian
    J_matrix = [x(2)-x(1), x(3)-x(1); y(2)-y(1), y(3)-y(1)];
    J = abs(det(J_matrix));
    JT_inv = inv(transpose(J_matrix));
    
    % Compute alok1
    alok1=zeros(3,3);
    
    for i = 1: 3
        for j = 1: 3
            I = 0;
            for n = 1: length(xi)
                % xn = x(1) + xi(n)*(x(2)-x(1)) + eta(n)*(x(3)-x(1));
                % yn = y(1) + xi(n)*(y(2)-y(1)) + eta(n)*(y(3)-y(1));
                
                temp = f(xi(n), eta(n));
                I = I + nu(n)*temp(i)*temp(j);
            end
            alok1(i, j) = I*J;
        end
    end
    
    % Compute alok2
    alok2 = zeros(3,3);
    for i = 1: 3
        for j = 1: 3
            I = 0;
            for n = 1: length(xi)
                % xn = x(1) + xi(n)*(x(2)-x(1)) + eta(n)*(x(3)-x(1));
                % yn = y(1) + xi(n)*(y(2)-y(1)) + eta(n)*(y(3)-y(1));
    
                I = I + nu(n)*dot((JT_inv*g(:, i)), JT_inv*g(:, j));
            end
            alok2(i, j) = I*J;
        end
    end