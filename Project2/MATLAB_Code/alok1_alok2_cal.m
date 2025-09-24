function [alok1, alok2] = alok1_alok2_cal(alok0, p1, p2, p3, q1, q2, q3)
    P = [p1, p2, p3];
    Q = [q1, q2, q3];

    Ap = 0.5*norm(cross(p2-p1, p3-p1));
    Aq = 0.5*norm(cross(q2-q1, q3-q1));

    alok1 = zeros(3);
    alok2 = zeros(3);

    for i = 1: 3
        for j = 1: 3
            p = i;
            p_prime = j;

            q = mod(i, 3)+1;
            q_prime = mod(j, 3)+1;
        
            ind = zeros(4);

            ind(1, 1) = q+2;
            ind(1, 2) = q+1;
            ind(1, 3) = q_prime+2;
            ind(1, 4) = q_prime+1;
            ind(2, 1) = q+2;
            ind(2, 2) = q+1;
            ind(2, 3) = p_prime+2;
            ind(2, 4) = p_prime+1;
            ind(3, 1) = p+2;
            ind(3, 2) = p+1;
            ind(3, 3) = q_prime+2;
            ind(3, 4) = q_prime+1;
            ind(4, 1) = p+2;
            ind(4, 2) = p+1;
            ind(4, 3) = p_prime+2;
            ind(4, 4) = p_prime+1;

            for m = 1: 4
                for n = 1: 4
                    if ind(m, n) > 3
                        ind(m, n) = ind(m, n)-3;
                    end
                end
            end
            
            alok1(i, j) = dot(P(:, ind(1, 1))-P(:, ind(1, 2)), Q(:, ind(1, 3))-Q(:, ind(1, 4)))*alok0(p, p_prime)...
                -dot(P(:, ind(2, 1))-P(:, ind(2, 2)), Q(:, ind(2, 3))-Q(:, ind(2, 4)))*alok0(p, q_prime)...
                -dot(P(:, ind(3, 1))-P(:, ind(3, 2)), Q(:, ind(3, 3))-Q(:, ind(3, 4)))*alok0(q, p_prime)...
                +dot(P(:, ind(4, 1))-P(:, ind(4, 2)), Q(:, ind(4, 3))-Q(:, ind(4, 4)))*alok0(q, q_prime);
            alok1(i, j) = alok1(i, j)/(4*Ap*Aq);
    
            alok2(i, j) = sum(alok0, "all")/(Ap*Aq);
        end
    end
end