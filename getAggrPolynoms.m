function F = getAggrPolynoms(x, AR_coeff, freq_TA, sampl_type)
    p = length(AR_coeff);    
    l = 0;

   if strcmp(sampl_type, 'flow') == 1
        A = zeros((1 + p) * freq_TA, 1 + p);
        for i = 1:(1 + p) * freq_TA
            for j = 1:1 + p
                if i == j
                    A(i, j) = 1;
                else if ((abs(i - j) <= (freq_TA - 1)*(p + 1)) && (i - j > 0))
                    A(i, j) = x(i - j);
                    l = i - j;
                    end
                end
            end
        end
        F1 = horzcat(1, AR_coeff);
        F1 = F1';
        F2 = zeros((1 + p) * freq_TA, 1);
        for i = 0:p
            for j = 1:freq_TA
                if i == 0
                    F2(i * freq_TA + j, 1) = 1;
                else
                    F2(i * freq_TA + j, 1) = x(l);
                end
            end
            l = l + 1;
        end
    else
        A = zeros(1 + p * freq_TA, 1 + p);
        for i = 1:1 + p * freq_TA
            for j = 1:1 + p
                if i == j
                    A(i, j) = 1;
                else if (abs(i - j) <= (p * (freq_TA - 1)) && i - j > 0)
                    A(i, j) = x(i - j);
                    l = i - j;
                    end
                end
            end
        end
        l = l + 1;
        F1 = horzcat(1, AR_coeff);
        F1 = F1';
        F2 = zeros(1 + p * freq_TA, 1);
        F2(1, 1) = 1;
        for i = 2:1 + p * freq_TA
            if mod(i - 1, freq_TA) == 0
                    F2(i, 1) = x(l);
                    l = l + 1;
            end
        end
    end
    F = A * F1 - F2;
end

