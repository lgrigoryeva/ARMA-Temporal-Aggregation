function F = getAggrPolynoms(x, AR_coeff, aggr_vector, K_ast)
    p = length(AR_coeff);
    
    K = length(aggr_vector(:));
    
    aggr_vector_refl = (flipud(aggr_vector(:)))';
    n = K * (p + 1) - K_ast - p;
    A = zeros(n + p + 1 + K_ast - 1, 1 + p);
    
    for i = 1:n + p + 1
        for j = 1:1 + p
            if i == j
                A(i, j) = x(1);
            else if ((abs(i - j) <= n) && (i - j > 0))
                    A(i, j) = x(i - j + 1);
                    
                end
            end
        end
    end
    l = n + 2;
    F1 = [1; AR_coeff(:)];
    
    F2 = zeros(K * (p + 1), 1);
    F2(1:K, 1) = aggr_vector_refl(:);
    
    for i = 2:p + 1
        F2((i - 1) * K + 1:i * K, 1) = aggr_vector_refl(:).*x(l);
        l = l + 1;
    end
    F = A * F1 - F2;
end

