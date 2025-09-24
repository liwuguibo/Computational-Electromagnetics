clc;
clear;

z = -10:0.05:10;
N =length(z);

epsilon_0 = 8.85e-12;

V1 = ones(1, N);

for i = 1 : N
    [V1(i), Q] = exercise14b_Tang_1(1, 2, z(i));
end

V_estimate = @(h) Q ./ (4*pi*epsilon_0*abs(h));

V2 = V_estimate(z);

plot(z, V1, 'b');
hold on;
plot(z, V2, 'r');
legend('Integral', 'Estimate');
xlabel('z/m');
ylabel('Potential/V');
hold off;

function [V, Q] = exercise14b_Tang_1(R1, R2, Z)

epsilon_0 = 8.85e-12;

% Define integral interval
a_rho = R1;
b_rho = R2;

% Generate points and weights
P = 4;
[xi, nu] = gausslegendre(P);

% Apply points and weights to variants
rho = a_rho + (b_rho-a_rho)*xi;
w_rho = (b_rho-a_rho)*nu;

f = @(x) 2.32.*(x.^4) - 13.75.*(x.^3) + 31.80.*(x.^2) -33.83.*x + 13.90;

g = f(rho).*rho./sqrt(Z^2 + rho.^2);

% Compute integral
V = sum(g.*w_rho);
V = V*2*pi/4/pi/epsilon_0/1e9;

Q = sum(f(rho).*w_rho.*rho);
Q = 2*pi*Q/1e9;
end