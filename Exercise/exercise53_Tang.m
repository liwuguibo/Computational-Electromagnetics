%% Clear worksparce and command window  
clc;
clear;
%% Defining parameters
gamma = 1.781072418;
a = 1/4;
k = pi/2;
Q = 3;

[xi, nu] = gausslegendre(Q);

x = a*xi;
w = a*nu;
    
%% Direct Gauss_Legendre
integrand = besselh(0, 1, k*x);
integral_gauss = sum(dot(w, integrand));
    
%% Sigularity subtraction
% Subtraction part
Hs = (2i/pi) * log(k*x*gamma/2);
corrected_integrand = integrand - Hs;
part1 = sum(dot(w, corrected_integrand));

% Analytical part
integral_log = a*(log(k*a)-1);
part2 = (2i/pi)*integral_log;

integral_singularity = part1+part2;
    
%% Output
disp(['Results with N=', num2str(Q)]);
disp(['Direct Gauss-Legendre Integral: ', num2str(integral_gauss)]);
disp(['Singularity Subtraction Integral: ', num2str(integral_singularity)]);