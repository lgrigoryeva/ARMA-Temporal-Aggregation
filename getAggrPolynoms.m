function F = getAggrPolynoms(x, AR, period_TA, aggr_type)
% It sets up Equation (2.7) in Paper 2
    p = length(AR);
    l = 0;

    if strcmp(aggr_type, 'flow') == 1
        A = zeros((1 + p) * period_TA, 1 + p);
        for i = 1:(1 + p) * period_TA
            for j = 1:1 + p
                if i == j
                    A(i, j) = 1;
                else if ((abs(i - j) <= (period_TA - 1)*(p + 1)) && (i - j > 0))
                        A(i, j) = x(i - j);
                        l = i - j;
                    end
                end
            end
        end
        F1 = horzcat(1, AR);
        F1 = F1';
        F2 = zeros((1 + p) * period_TA, 1);
        for i = 0:p
            for j = 1:period_TA
                if i == 0
                    F2(i * period_TA + j, 1) = 1;
                else
                    F2(i * period_TA + j, 1) = x(l);
                end
            end
            l = l + 1;
        end
    else
        A = zeros(1 + p * period_TA, 1 + p);
        for i = 1:1 + p * period_TA
            for j = 1:1 + p
                if i == j
                    A(i, j) = 1;
                else if (abs(i - j) <= (p * (period_TA - 1)) && i - j > 0)
                        A(i, j) = x(i - j);
                        l = i - j;
                    end
                end
            end
        end
        l = l + 1;
        F1 = horzcat(1, AR);
        F1 = F1';
        F2 = zeros(1 + p * period_TA, 1);
        F2(1, 1) = 1;
        for i = 2:1 + p * period_TA
            if mod(i - 1, period_TA) == 0
                F2(i, 1) = x(l);
                l = l + 1;
            end
        end
    end
    
    % constructs the target function F which has to satisfy F = 0
    F = A * F1 - F2;
end

