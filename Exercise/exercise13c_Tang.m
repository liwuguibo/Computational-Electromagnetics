clc;
clear;

% Define integral interval
a_x = 0; 
b_x = 3*pi/2;

% Generate points and weights
P = 19;
[xi, nu_x] = gausslegendre(P);

% Apply points and weights to variant
x = a_x + (b_x-a_x) * xi;
w_x = (b_x-a_x) * nu_x;

f = 3*(sin(2*x).^2) - 2*cos(1.5*(x.^2));

% Compute integral
I = sum(f.*w_x);

disp(['Results: ', num2str(I)]);
disp([num2str(P) ,' integration points are needed']);