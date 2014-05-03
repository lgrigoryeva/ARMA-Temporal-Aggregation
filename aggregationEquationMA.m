function [ f ] = aggregationEquationMA(x, q_new, freq_TA, coefficients_C_L, coefficients_dC_L, MA_TA_L, var_e, var_e_TA)

    %introduce shift matrix
    autocov_vector_initial = zeros(q_new + 1, 1);
    autocov_vector_TA = zeros(q_new + 1, 1);
    [~, n] = size(coefficients_C_L);
    coefficients_dMA_B(1) = 0;
        k = 2;
        for i = 1:q_new
            coefficients_dMA_B(k) = x(i);
            k = k + 1;
        end
        coefficients_MA_B = cell2mat(toCellArray(MA_TA_L)); 
        mm = length(coefficients_MA_B);
        if mm < q_new + 1
            for i = 1:q_new + 1 - mm
                coefficients_MA_B = [coefficients_MA_B 0];
            end
        end
    for i = 1:q_new + 1
        S_k = getShiftMatrix(n, (i - 1) * freq_TA);
        autocov_vector_initial(i) = var_e * (coefficients_dC_L * S_k * coefficients_C_L' + coefficients_C_L * S_k * coefficients_dC_L');
        
        
        S_i = getShiftMatrix(q_new + 1, i - 1);
        
        
        autocov_vector_TA(i) = x(q_new + 1) * (coefficients_MA_B * S_i * coefficients_MA_B') + var_e_TA * (coefficients_dMA_B * S_i * coefficients_MA_B' + coefficients_MA_B * S_i * coefficients_dMA_B');
    end
    
    
    f1 = autocov_vector_TA;
    f2 = autocov_vector_initial;
    f = (f1 - f2)';
    
end

