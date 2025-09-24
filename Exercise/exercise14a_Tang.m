clc;
clear;

result = exercise14a_Tang_1(2, 0.25, 4);

disp(['The length of spiral is ', num2str(result)]);

function arclen = exercise14a_Tang_1(a, k, P)
    interval_l = 0;
    interval_u = 2*pi;

    x = @(t) a.*exp(k.*t).*cos(t);
    y = @(t) a.*exp(k.*t).*sin(t);
    
    h = 0.001;

    [xi, nu] = gausslegendre(P);
    
    phi = interval_l + (interval_u-interval_l) * xi;
    w = (interval_u-interval_l) * nu;

    g_x = @(t) (x(t+h) - x(t-h)) / (2*h);
    g_y = @(t) (y(t+h) - y(t-h)) / (2*h);

    g = sqrt(g_x(phi).^2 + g_y(phi).^2);

    arclen = sum(g.*w);
   
    phi_plot = 0: 0.01: 2*pi;
    plot(x(phi_plot), y(phi_plot));
    hold on;

    x_scatter = x(phi);
    y_scatter = y(phi);

    scatter(x_scatter, y_scatter, 'red', 'filled');
    title('The Curve of logarithmic spiral');
    hold off;
end