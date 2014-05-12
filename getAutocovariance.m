function autocov_vector = getAutocovariance(lag, MA, period_TA)
% Computes the autocovariance of the given lag and aggregation period_TA
% used for example in the right hand side of Equation (2.21) in Paper 2
    autocov_vector = [];
    for i = 0:lag
        s = 0;
        for j = 1:length(MA) - i * period_TA
            s = s + MA(j) * MA(j + i * period_TA);
        end
        autocov_vector = [autocov_vector s]; 
    end
end

