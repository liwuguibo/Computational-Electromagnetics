function [Exy, Hz] = exercise44_Tang(xy, p1, p2, p3, c, s, freq, mur)

    c_new = [c(2) c(3) c(1)];
    s_new = [s(2) s(3) s(1)];

    x_p = [p1(1) p2(1) p3(1)];
    y_p = [p1(2) p2(2) p3(2)];

    J_matrix = [x_p(2)-x_p(1), x_p(3)-x_p(1);...
        y_p(2)-y_p(1), y_p(3)-y_p(1)];
    J = abs(det(J_matrix));
    % JT_inv = inv(transpose(J_matrix));


    gradient_lambda1 = [(y_p(2) - y_p(3)), (x_p(3) - x_p(2))] ;
    gradient_lambda2 = [(y_p(3) - y_p(1)), (x_p(1) - x_p(3))] ;
    gradient_lambda3 = [(y_p(1) - y_p(2)), (x_p(2) - x_p(1))] ;


    lambda1 = ((y_p(2)-y_p(3)) * (xy(1)-x_p(3)) + (x_p(3)-x_p(2)) * (xy(2)-y_p(3))) / J;
    lambda2 = ((y_p(3)-y_p(1)) * (xy(1)-x_p(3)) + (x_p(1)-x_p(3)) * (xy(2)-y_p(3))) / J;
    lambda3 = 1 - lambda1 - lambda2;

    N1 = lambda2*gradient_lambda3 - lambda3*gradient_lambda2;
    N2 = lambda3*gradient_lambda1 - lambda1*gradient_lambda3;
    N3 = lambda1*gradient_lambda2 - lambda2*gradient_lambda1;


    E = s_new(1) * c_new(1) * N1 + s_new(2) * c_new(2) * N2 + s_new(3) * c_new(3) * N3;
    Exy = E.';

    omega = 2*pi*freq*1e9;
    mu0 = 4*pi*1e-7;
    mu = mur*mu0;
    curlE = 2*J^(-1)*(s_new(1)*c_new(1) + s_new(2)*c_new(2) + s_new(3)*c_new(3));
    Hz = -curlE / (1i * omega * mu);
end