clc;
clear;

% Define integral interval
a_x = -sqrt(2)/2; 
b_x = sqrt(2)/2;


% Generate points and weights
P = 6;
[xi, nu_x] = gausslegendre(P);

% Apply points and weights to variant
x = a_x + (b_x-a_x) * xi;
w_x = (b_x-a_x) * nu_x;

f = @(x) sqrt(1-x.^2);
h = 0.0001;

g = sqrt(1 + ((f(x+h) - f(x-h)) / (2*h)).^2);

% Compute integral
I = sum(g.*w_x);

disp(['Results: ', num2str(I)]);
disp([num2str(P) ,' integration points are needed']);