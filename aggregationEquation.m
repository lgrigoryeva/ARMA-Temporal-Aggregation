function [ f ] = aggregationEquation(x, q_new, coefficients_TA,  var_e, freq_TA)
    
    aggr_polynom = 1;
    
    for i = 2:q_new + 1
       aggr_polynom(i) = x(i - 1); 
    end
    
    autocov_vector_TA = getAutocovariance(q_new, aggr_polynom,1);
    autocov_vector_TA = autocov_vector_TA*x(q_new + 1);
    
    autocov_vector_initial = getAutocovariance(q_new, coefficients_TA,freq_TA);
    autocov_vector_initial = autocov_vector_initial*var_e;
    
    f1 = autocov_vector_TA;
    f2 = autocov_vector_initial;
    f = (f1 - f2)';
    
end

