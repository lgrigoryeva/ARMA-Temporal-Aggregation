function S_k = getShiftMatrix(n, shift)
% This function produces the lower-shift matrix S_k of dimension n x n
% with the given shift value (see page 11 of Paper 2)

    S_k = zeros(n);
    
    for i = 1:n
        for j = 1:n
            if (i - j) == shift
                S_k(i, j) = 1;
            end
        end
    end
    
