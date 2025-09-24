function I = exercise23a_Tang(p1, p2, p3, Q)

% Define integral function
f = @(x, y) 2*(x.^2)/3 + 3*y/2;

% Compute integration points and weights
[xi, eta, nu] = gausstriangle(Q);

% Compute Jacobian
J = abs((p2(1)-p1(1))*(p3(2)-p1(2)) - (p3(1)-p1(1))*(p2(2)-p1(2)));

% Compute numerical integral
I = 0;

for i = 1: length(xi)
    x = p1(1) + xi(i)*(p2(1)-p1(1)) + eta(i)*(p3(1)-p1(1));
    y = p1(2) + xi(i)*(p2(2)-p1(2)) + eta(i)*(p3(2)-p1(2));

    I = I + nu(i)*f(x,y);
end

I = I*J;
