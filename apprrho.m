function appr_rho = apprrho(A,B)
initial_state = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
all_operator = B * A;
q2rho = (all_operator) * initial_state * (all_operator)';
appr_rho = zeros(4,4);
I = eye(2);
X = [0 1; 1 0];
Y = [0 -1i; 1i 0];
Z = [1 0; 0 -1];

for i = 1:4
    if i == 1
        lambda = I;
    end
    if i == 2
        lambda = X;
    end
    if i == 3
        lambda = Y;
    end
    if i == 4
        lambda = Z;
    end
    appr_rho = appr_rho + kron(I, lambda) * trace(kron(I, lambda) * q2rho);
end
end

        
