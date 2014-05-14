function [AR_H, MA_H, sigma2_H] = temporalAggregationForARMA(p, q, Coeff, w, structInitCond)
    %
    %PARAMETERS INITIALIZATION FOR TEMPORARY AGGREGATION PROCEDURE
    %
    AR_coeff = - Coeff.AR;
    
    if isfield(Coeff, 'MA') && q
        MA_coeff = Coeff.MA;
    else
        MA_coeff = 0;
    end
    
    sigma2 = Coeff.K;
    
    K = length(w);
    
    K_ast = 1;
    while w(K_ast) == 0
        K_ast = K_ast + 1;
    end
    
    n = K * (p + 1) - K_ast - p;
    
    x_length = n + p + 1;
    
    %
    %TEMPORAL AGGREGATION PROCEDURE
    %        
    x0 = -1 * ones(x_length, 1);
      
    %get the T(L) aggregation lag operator and A*(L)=T(L)*A(L) lag operator for aggregated AR
    %create auxiliary function
    f = @(x)getAggrPolynoms(x, AR_coeff, w, K_ast);
    %set up options for the solver
    options = optimset('Display', 'iter', 'Algorithm', 'trust-region-reflective', 'TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 2000000, 'MaxIter', 2000000);
    %get the solution for the coefficients
    [x, ~, exitflag, output] = fsolve(f, x0, options);
    if exitflag <= 0
        display('problem when solving getAggrPolynoms');
        display(output);
        return;
    end
    
    %
    %PROCESSING THE RESULT x
    %    
    T_L = LagOp({1});
    coef_TA(1) = 1;
    AR_H = 0;
    
    k = 1;
    for i = 2:n + 1
        T_L.Coefficients{k} = x(i);
        coef_TA(k + 1) = x(i);
        k = k + 1;
    end
    
    k = 2;
    for i = n + 2:length(x)
        AR_H(k - 1) = - x(i);
        k = k + 1;
    end
    
    %new order of MA of temporarily aggregated ARMA
    q_new = idivide(int16(K * (p + 1) + q - p - K_ast), K);

    %get the new MA aggregated polynomial
    if q
        MA_L = mtimes(T_L, LagOp([1 MA_coeff]));
        [coefficients_H] = toCellArray(MA_L);
    else
        [coefficients_H] = toCellArray(T_L);
    end

    %
    %TEMPORAL AGGREGATION LOWERING THE ORDER OF L^k = B
    %
    %define the variance of the TA model
    coefficients_H = cell2mat(coefficients_H);
    %set up initial point for the solver
    
    if nargin < 5
        x0 = zeros(q_new + 1, 1);
    else
        x0 = [structInitCond.MA, structInitCond.K];
    end
    %set the auxiliary function for estimation of the parameters of the new
    %model
    f = @(y)aggregationEquation(y, q_new, coefficients_H, sigma2, K);
    %set up the options for the solver
    options = optimset('Display', 'iter', 'Algorithm','trust-region-dogleg', 'TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 500000, 'MaxIter', 200000);
    [y, ~, exitflag, output] = fsolve(f, x0, options);
    if exitflag <= 0
        display('problem when solving aggregationEquation');
        display(output);
        x0 = x0*0;
        [y, ~, exitflag, output] = fsolve(f, x0, options);
        if exitflag <= 0
            display('problem when solving aggregationEquation');
            display(output);
            return;
        end
    end

    %y contains parameters of the aggregated model and variance of the
    %aggregated residuals

    MA_H = [0];
    k = 1;
    for i = 1:q_new
        MA_H(k) = y(i);
        k = k + 1;
    end

    sigma2_H = y(q_new + 1);

end

    

