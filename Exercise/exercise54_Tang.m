function sigma = exercise54_Tang(x, coord, elements)
    % Constants
    c0 = physconst('LightSpeed');
    f = 1e8;
    lambda = c0/f;
    k = 2*pi/lambda;
    r_obs = 100*lambda;

    % Observation angles
    N_phi = 360;
    phi = linspace(0, 2*pi, N_phi);
    x_obs = r_obs * cos(phi);
    y_obs = r_obs * sin(phi);

    N_elements = length(elements(1, :));

    elem_centers = zeros(2, N_elements);
    elem_lengths = zeros(1, N_elements);

    for n = 1: N_elements
        p1 = coord(:, elements(1, n));
        p2 = coord(:, elements(2, n));
        
        elem_centers(:,n) = (p1 + p2)/2;
        elem_lengths(n) = norm(p2 - p1);
    end

    % Compute scattered field at each observation point
    Ezsca = zeros(1, N_phi);
    for m = 1: N_elements
        obs_point = [x_obs(m); y_obs(m)];

        for n = 1: N_elements
            r_prime = elem_centers(:, n);
            R = norm(obs_point - r_prime);
            G = 1i/4*besselh(0, 1, k*R);
            Ezsca(m) = Ezsca(m) + G*x(n)*elem_lengths(n);
        end
    end

    % Compute bistatic RCS
    sigma = 2*pi*r_obs*(abs(Ezsca).^2);

    % Normalize and convert to dB scale
    sigma_norm = sigma / max(sigma);
    sigma_dB = 10 * log10(sigma_norm);

    % Plot
    figure;
    plot(rad2deg(phi), sigma_dB);
    xlabel('\phi (degrees)');
    ylabel('Normalized RCS (dB)');
    title('Normalized Bistatic RCS in dB');
    grid on;
end