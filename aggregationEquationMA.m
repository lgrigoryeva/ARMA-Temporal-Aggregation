function [ f ] = aggregationEquationMA(x, q_star, period_TA, C_L, dC_L, Theta_star_B_polynom, sigma2, sigma2_TA)

% This function sets up the equation (2.27) in Paper 2


    % C_L = \bar{T(L) * Theta(L)} in Equation (2.27) of Paper 2
    % dC_L = \bar{dT(L)/dbeta_i * Theta(L)} + \bar{T(L) * dTheta(L)/dbeta_i}  in Equation (2.27) of Paper 2
    % Theta_star_B_polynom - Theta_star_B in LagOp form
    
    autocov_vector_initial = zeros(q_star + 1, 1);
    autocov_vector_TA = zeros(q_star + 1, 1);
    
    [~, n] = size(C_L);
    
    dTheta_star_B = zeros(q_star + 1, 1);
    dTheta_star_B(1) = 0;
    
    k = 2;
    for i = 1:q_star
        dTheta_star_B(k) = x(i);
        k = k + 1;
    end
    
    Theta_star_B = cell2mat(toCellArray(Theta_star_B_polynom));
    mm = length(Theta_star_B);
    if mm < q_star + 1
        for i = 1:q_star + 1 - mm
            Theta_star_B = [Theta_star_B 0];
        end
    end
    
    % construct the equations
    for i = 1:q_star + 1
        %introduce shift matrix
        % left hand side of Equation (2.27) of Paper 2
        S_k = getShiftMatrix(n, (i - 1) * period_TA);
        autocov_vector_initial(i) = sigma2 * (dC_L * S_k * C_L' + C_L * S_k * dC_L');
        % right hand side of Equation (2.27) of Paper 2
        S_i = getShiftMatrix(q_star + 1, i - 1);
        autocov_vector_TA(i) = x(q_star + 1) * (Theta_star_B * S_i * Theta_star_B') + sigma2_TA * (dTheta_star_B * S_i * Theta_star_B' + Theta_star_B * S_i * dTheta_star_B');
    end
    
    % constructs the target function f which has to satisfy f = 0
    f1 = autocov_vector_TA;
    f2 = autocov_vector_initial;
    f = (f1 - f2)';
    
end

