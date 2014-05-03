function S_k = getShiftMatrix(n, shift)

    S_k = zeros(n);
    
    for i = 1:n
        for j = 1:n
            if (i - j) == shift
                S_k(i, j) = 1;
            end
        end
    end
    
