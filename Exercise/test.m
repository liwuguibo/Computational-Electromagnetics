%% EX33
clc;
clear;

load("exercise33_mesh.mat");

U = 5;
epsr = 2.3;
x = exercise33_Tang(coord, elements, nodes, U, epsr);
plot_scalar2D(x, coord, elements);

%% EX34
clc;
clear;

load("exercise34_mesh.mat");
[fc, xc] = exercise34_Tang(coord, elements, nodes);
plot_scalar2D(xc(:, 6), coord*1000, elements);

%% EX43
clc;
clear;

p1 = [-0.2; -0.3];
p2 = [1.5; 0.5];
p3 = [1.0; 1.5];

Q = 3;
[alok3, alok4] = exercise43_Tang(p1, p2, p3 ,Q);

%% EX44
clc;
clear;

[xi, eta, nu] = gausstriangle(10);

p1 = [0; 0];
p2 = [2; 0];
p3 = [0; 1];
c = [0, 0, 1];
s = c;

xy = p1 + (p2 - p1)*xi + (p3 - p1)*eta;

N = length(xy);

Exy = zeros(2, N);
Hz = zeros(2, N);

for i = 1: N
    [Exy(:, i), Hz(:, i)] = exercise44_Tang(xy(:, i), p1, p2, p3, c, s, 1, 1);
end

X = xy(1, :);
Y = xy(2, :);
U = Exy(1, :);
V = Exy(2, :);

figure;
quiver(X, Y, U, V);
hold on;
fill([p1(1), p2(1), p3(1)], [p1(2), p2(2), p3(2)], 'k', 'FaceAlpha', 0.1);

%% alok6
clc;
clear;

load("project_work1_mesh.mat");

N = 1;
vertices = elements(:, N);
p1 = coord(:, vertices(1));
p2 = coord(:, vertices(2));
p3 = coord(:, vertices(3));

Q = 3;

alok6 = project_work11_alok6_Tang(p1, p2, p3, Q);

%% Exercise 54
clc;
clear;

load("exercise54_data.mat");
sigma = exercise54_Tang(x, coord, elements);

%% Exercise 63 and 64
clc;
clear;

load("exercise63_mesh.mat");
x = exercise63_Tang(coord, elements);

theta = pi/2;
z = 0;
N_field = 100;

a = 0.81e-3;
b = 2.9e-3;
r = linspace(a, b, N_field);

xyz = zeros(3, N_field);

for i = 1: N_field
    xyz(:,i) = [r(i)*cos(theta); r(i)*sin(theta); z];
end

phi = exercise64_Tang(xyz, coord, elements, x);

figure;
plot(r, phi);
xlabel('Radius/m');
ylabel('Potential/V');
title('Potential v.s radius in a coaxial cable');

%% Exercise 73
clc;
clear;

load('exercise73_data.mat');
P = 3;
freq = 1e9;
c0 = physconst('LightSpeed');
lambda = c0/freq;
L = lambda/2;

F = exercise73_Tang(freq, L, P, XYZ);

plot_radiation_pattern(F, XYZ, ET);

%% Exercise 83
clc; 
clear;

k = 0.7;

p1 = [0.7; 0.5; 0.1];
p2 = [0.4; 0.7; 0.1];
p3 = [0.2; 0.5; 0.9];

q1 = [0.8; 1.9; 0.9];
q2 = [1.5; 1.6; 0.8];
q3 = [1.6; 1.7; 1.9];

Q = 3;

[xi, eta, nu] = gausstriangle(Q);

alok0 = exercise83_Tang(k, p1, p2, p3, q1, q2, q3, xi, eta, nu);

disp(alok0);

%% Exercise 84
clc; 
clear;

k = 0.25;

p1 = [0.8; 0.9; 0.1];
p2 = [0.9; 0.6; 0.1];
p3 = [0.3; 0.4; 0.9];

q1 = [0.8; 0.9; 0.1];
q2 = [0.9; 0.6; 0.1];
q3 = [0.9; 0.5; 0.75];

Q = 5;

[xi, eta, nu] = gausstriangle(Q);

alok0 = exercise84_Tang(k, p1, p2, p3, q1, q2, q3, xi, eta, nu);

disp(alok0);