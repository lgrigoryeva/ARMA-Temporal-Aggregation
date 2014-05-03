function [V_beta, gamma_array_p, gamma_array_q, gamma_array_pq] = getAsymptMLEcovar(AR, MA, precision, MA_el)
    
    p = length(AR);
    q = length(MA);
    V_beta = zeros(1);
    
    if p
    
        gamma_vector_p = zeros(p, 1);
        for i = 1:p
            gamma_vector_p(i) = getAutocovar(AR, [], i - 1, precision, MA_el);
        end

        gamma_array_p = zeros(p);
        for i = 1:p
            for j = 1:p
                gamma_array_p(i, j) = gamma_vector_p(abs(i - j) + 1);
            end
        end
        if ~q
            V_beta = inv(gamma_array_p);
            return
        end
    end
    
    if q
        
        gamma_vector_q = zeros(q, 1);
        for i = 1:q
            gamma_vector_q(i) = getAutocovar(-MA, [], i - 1, precision, MA_el);
        end

        gamma_array_q = zeros(q);
        for i = 1:q
            for j = 1:q
                gamma_array_q(i, j) = gamma_vector_q(abs(i - j) + 1);
            end
        end
        if ~p
            V_beta = inv(gamma_array_q);
            return
        end
    end
    
    if p && q
        gamma_array_pq = getAutocovarMixed(AR, -MA, precision, MA_el);

        V_beta = inv([gamma_array_p, gamma_array_pq; gamma_array_pq', gamma_array_q]);
    end
end