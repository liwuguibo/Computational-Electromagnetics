clc;
clear;

% Define integral interval
a = -1; 
b = 3;

% Generate points and weights
[xi, nu] = gausslegendre(3);

% Apply points and weights to variant
x = a + (b-a)*xi;
w = (b-a)*nu;
f = (x.^5)/3 - 2*(x.^3) + x - 1;

% Compute integral
I = sum(f.*w);

disp(['Results: ', num2str(I)]);
disp('3 integration points are needed');