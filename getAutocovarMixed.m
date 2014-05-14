function [ gamma_array ] = getAutocovarMixed(AR, MA, precision, recursStartLevel)
% This function is used in the function getAsymptMLEcovar for computing the
% (1,2),(2,1) blocks of SIGMA_beta, see page 7 of Paper 2
    p = length(AR);
    q = length(MA);
    [Psi_p, Psi_q] = getAutocovMixedRecur(AR, MA, max(p, q), precision, recursStartLevel);
    
    gamma_array = zeros(p, q);
    for i = 1:p
        for j = 1:q
            if (i < j)
                gamma_array(i, j) = sum(reshape(Psi_p((j - i + 1):(end)), [], 1).*reshape(Psi_q(1:(end - j + i)), [], 1), 'double');
            else
                gamma_array(i, j) = sum(reshape(Psi_p(1:(end - i + j)), [], 1).*reshape(Psi_q((i - j + 1):end), [], 1), 'double');
            end
        end
    end
    
end

