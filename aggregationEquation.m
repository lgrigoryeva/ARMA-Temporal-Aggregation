function [ f ] = aggregationEquation(x, q_star, params_TA, sigma2, period_TA)
% Sets up Equation (2.22) in Paper 2
    
    aggr_polynom = zeros(q_star + 1, 1);
    aggr_polynom(1) = 1;
    
    for i = 2:q_star + 1
       aggr_polynom(i) = x(i - 1); 
    end
    
    % right hand side of Equation (2.22)
    autocov_vector_TA = getAutocovariance(q_star, aggr_polynom, 1);
    autocov_vector_TA = autocov_vector_TA * x(q_star + 1);
    
    % left hand side of Equation (2.22)
    autocov_vector_initial = getAutocovariance(q_star, params_TA, period_TA);
    autocov_vector_initial = autocov_vector_initial * sigma2;
    
    f1 = autocov_vector_TA;
    f2 = autocov_vector_initial;
    
    % constructs the target function f which has to satisfy f = 0
    f = (f1 - f2)';
    
end

