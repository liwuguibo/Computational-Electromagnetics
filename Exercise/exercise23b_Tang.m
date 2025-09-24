function Qtot = exercise23b_Tang(elements, coord, c, Q)

    epsilon0 = 8.854e-12;
    coeff = 1/(4*pi*epsilon0);
    
    I = 0;
    
    for i = 1: size(elements, 2)
        % Get the nodes
        nodes = elements(:,i);
        
        % Get the positions of the nodes
        x = coord(1, nodes);
        y = coord(2, nodes);
        
        % Get local c value
        c_local = c(nodes);
        
        % Compute integration points and weights
        [xi, eta, nu] = gausstriangle(Q);
        
        % Define integral function
        f = @(x, y) c_local(1)*(1-x-y) + c_local(2)*x + c_local(3)*y;
        
        % Compute Jacobian
        J = abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)));
        
        % Compute integral
        Ii = 0;
    
        for n = 1 : length(xi)
            xn = x(1) + xi(n)*(x(2)-x(1)) + eta(n)*(x(3)-x(1));
            yn = y(1) + xi(n)*(y(2)-y(1)) + eta(n)*(y(3)-y(1));
    
            Ii = Ii + nu(n)*f(xn, yn);
        end
        Ii = Ii*J;
        I = I + Ii;
    end
    
    Qtot = coeff*I;