clc;
clear;

N_ports = 3;

load('S_dB_Nf=101.mat');

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