clc;
clear;

% Define integral interval
a_x = 0; 
b_x = 4;

a_y = -1;
b_y = 3;

% Generate points and weights
P = 4;
[xi, nu_x] = gausslegendre(P);
[yi, nu_y] = gausslegendre(P);

% Apply points and weights to variant
x = a_x + (b_x-a_x) * xi;
w_x = (b_x-a_x) * nu_x;

y = a_y + (b_y - a_y) * yi;
w_y = (b_y - a_y) * nu_y;

[X,Y] = meshgrid(x,y);

f = X.*Y.^2 - 3*X.^3 + X.*Y.^7 -3;

% Compute integral
I = w_y*(f*w_x');

disp(['Results: ', num2str(I)]);
disp([num2str(P^2) ,' integration points are needed']);