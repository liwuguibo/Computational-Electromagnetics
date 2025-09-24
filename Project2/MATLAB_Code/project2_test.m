%% alok1 alok2 test T1 T1 singular
clc;
clear;
load('project_work2_data.mat');
load('project_work2_mesh.mat');

P2 = 4;

freq = 1e9;

c0 = physconst('LightSpeed');
lambda = c0/freq;
k_wave = 2*pi/lambda;

node_indices_1 = elements(:, 1);
node_indices_2 = elements(:, 1);

p1 = coord(:, node_indices_1(1));
p2 = coord(:, node_indices_1(2));
p3 = coord(:, node_indices_1(3));

q1 = coord(:, node_indices_2(1));
q2 = coord(:, node_indices_2(2));
q3 = coord(:, node_indices_2(3));

alok0 = alok0_cal_singular(k_wave, p1, p2, p3, P2);

[alok1, alok2] = alok1_alok2_cal(alok0, p1, p2, p3, q1, q2, q3);

format long e;
disp(alok1);
disp(alok2);

%% alok1 alok2 test T1 T2 near-singular
node_indices_1 = elements(:, 1);
node_indices_2 = elements(:, 2);

p1 = coord(:, node_indices_1(1));
p2 = coord(:, node_indices_1(2));
p3 = coord(:, node_indices_1(3));

q1 = coord(:, node_indices_2(1));
q2 = coord(:, node_indices_2(2));
q3 = coord(:, node_indices_2(3));

alok0 = alok0_cal_near_singular(k_wave, p1, p2, p3, q1, q2, q3, P2);

[alok1, alok2] = alok1_alok2_cal(alok0, p1, p2, p3, q1, q2, q3);

clc;
format long e;
disp(alok1);
disp(alok2);

%% Project21
clc;
clear;

load('project_work2_mesh.mat');

P1 = 3;
P2 = 4;
freq = 1e9;

A = project_work21_Tang(P1, P2, coord, elements, basis, freq);

clc;
format long e;
disp(A(1:4, 1:2));
disp(A(8:11, 1:2));

%% Project22
clc;
clear;

load('project_work2_mesh.mat');

f_min = 0.5e9;
f_max = 1.5e9;
Nf = 101;
freq_list = linspace(f_min, f_max, Nf);

N_ports = length(ports);
S_dB_all = zeros(N_ports, N_ports, Nf);

for idx = 1: Nf
    freq = freq_list(idx);
    S = project_work22_Tang(coord, elements, basis, ports, freq);
    S_dB_all(:, :, idx) = 20*log10(abs(S));
end

close all;

figure(1);
for i = 1: N_ports
    for j = 1: N_ports
        k = (i-1)*N_ports+j;
        subplot(N_ports, N_ports, k);
        plot(freq_list/1e9, squeeze(S_dB_all(i, j,:)), 'LineWidth', 1.2);
        xlabel('Frequency (GHz)');
        ylabel('Magnitude (dB)');
        title(['S_{' num2str(i) num2str(j) '}']);
        grid on;
        xlim([f_min f_max]/1e9);
    end
end
sgtitle('3×3 S-Parameters');

%% dlok T1
clc;
clear;

load('project_work2_mesh.mat');
load("project_work2_data.mat");

p1 = coord(:, elements(1, 1));
p2 = coord(:, elements(2, 1));
p3 = coord(:, elements(3, 1));

r = XYZ(:, 1);

P = 3;

freq = 1e9;
c0 = physconst('LightSpeed');
lambda = c0/freq;
k_wave = 2*pi/lambda;

dlok = dlok_cal(P, p1, p2, p3, r, k_wave);

clc;
format long e;
disp(dlok);

%% Project 23
clc;
clear;

load('project_work2_mesh.mat');
load("project_work2_data.mat");

P1 = 3;
P2 = 4;
P = 3;
freq = 0.96e9;

A = project_work21_Tang(P1, P2, coord, elements, basis, freq);

N_interior = max(max(abs(basis)));
N_ports = length(ports);

B = zeros(N_interior, 1);
phaseshift = deg2rad(0);

for m = 1: N_ports
    B(ports(m), 1) = exp(1i*phaseshift*(m-1));
end

x = A\B;

F = project_work23_Tang(P, XYZ, coord, elements, basis, freq, x);

close all;

plot_radiation_pattern(F, XYZ, ET);
