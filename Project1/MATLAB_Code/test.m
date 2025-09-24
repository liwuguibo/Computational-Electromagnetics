%% EX23 function test
clc;
clear;

p1 = [0.2; 0.8];
p2 = [1; 0.5];
p3 = [0.5; 1.5];

Q = 3;

[alok1, alok2] = exercise24_Tang(p1, p2, p3, Q);

% disp(alok1);
% disp(alok2);

%% EX43 function test
clc;
clear;

p1 = [-0.2; -0.3];
p2 = [1.5; 0.5];
p3 = [1.0; 1.5];

Q = 3;
[alok3, alok4] = exercise43_Tang(p1, p2, p3 ,Q);

% disp(alok3);
% disp(alok4);

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

[alok1, alok2, alok3, alok4, alok6] = calc_aloks_Tang(p1, p2, p3, Q);

% disp(alok6);

%% 1
clc;
clear;

load("project_work1_mesh.mat");
f0 = 5e10;

[A, B, C, D] = project_work11_Tang(coord, elements, basis, edges, nodes, f0);

% disp(A(1:4, 1:4));
% disp(B(1:4, 1:4));
% disp(C(1:4, 373:376));
% disp(D(1:4, 1:4));

%% 2
clc;

[beta, fc, xET, xEZ] = project_work12_Tang(A, B, C, D, edges, nodes, f0);

%% 3
clc;

N = 4;
xETj = xET(:, N);
xEzj = xEZ(:, N);
betaj = beta(N);

[E, H, xy] = project_work13_Tang(xETj, xEzj, coord, elements, basis, f0, betaj);

%% ET Plot

plot_constant2D(E(1:2,:), coord, elements);

%% HT Plot

plot_constant2D(H(1:2,:), coord, elements);

%% E_Z Plot

plot_scalar2D(E(3,:), coord, elements);

%% H_Z Plot

plot_scalar2D(H(3,:), coord, elements);

%% Quiver Plot

figure(1);
quiver(xy(1,:), xy(2,:), real(E(1,:)), real(E(2,:)));
axis equal; 
title('Real part of Transverse Electric Field');

% figure(2);
% quiver(xy(1,:), xy(2,:), real(H(1,:)), real(H(2,:)));
% axis equal; 
% title('Real part of Transverse Magnetic Field');

%% Transverse E Plot All
clc;

figure(1);
sgtitle('Real part of Transverse Electric Field');  

for i = 1: 6
    xETj = xET(:, i);
    xEzj = xEZ(:, i);
    betaj = beta(i);

    [E, H, xy] = project_work13_Tang(xETj, xEzj, coord, elements, basis, f0, betaj);

    subplot(3,2,i);
    quiver(xy(1,:), xy(2,:), real(E(1,:)), real(E(2,:)));
    axis equal; 
    
    switch i
        case 1
            title('TE10');
        case 2
            title('TE20');
        case 3
            title('TE01');
        case 4
            title('TE11');
        case 5
            title('TM11');
        otherwise
            title('TE30');
    end
end