function [AR_H, MA_H, sigma2_H, T_L, dBeta] = temporalAggregation(p, q, Coeff, w, structInitCond)
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
    
    period_TA = length(w);
    K = period_TA;

    K_ast = 1;
    while w(K_ast) == 0
        K_ast = K_ast + 1;
    end
    
    n = K * (p + 1) - K_ast - p;
    
    x_length = n + p + 1;
    
    if K_ast == period_TA
        aggrType = 'stock';
    elseif K_ast == 1
        aggrType = 'flow';
    end
    
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
    
    %new order of MA part of temporarily aggregated ARMA
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
    
    if period_TA == 1
            x0 = [MA_coeff, sigma2];
        else
            x0 = [structInitCond.MA, structInitCond.K];
            if sum(x0) == 0
                x0 = [0 * structInitCond.MA, sigma2 * period_TA];
            end
    end
    %set the auxiliary function for estimation of the MA coefficients of the new
    %model
    f = @(y)aggregationEquation(y, q_new, coefficients_H, sigma2, period_TA);
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
    % if we need to compute the Jacobians
    if nargout > 4
        if strcmp(aggrType, 'flow')
            x0 = -10 * ones(((p + 1) * (period_TA - 1)) + p, 1);
        else
            x0 = -10 * ones((p * (period_TA - 1)) + p, 1);
        end
        derPolynom = [];
        for index = 1:p

            f = @(x)getAggrPolynomsDeriv(x, AR_coeff, coef_TA, period_TA, aggrType, index);
            options = optimset('Display', 'iter', 'Algorithm', 'trust-region-reflective', 'TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 200000, 'MaxIter', 200000);
            [x, ~, exitflag, output] = fsolve(f, x0, options);
            
            if exitflag <= 0
                display('problem when solving getAggrPolynomsDeriv');
                display(output);
                return;
            end
            
            
            derPolynom(index).dT_L = LagOp({0});
            derPolynom(index).dAR_pol = [0];
            k = 1;
            if strcmp(aggrType, 'flow')
                for i = 1:(p + 1) * (period_TA - 1)
                    derPolynom(index).dT_L.Coefficients{k} = x(i);
                    k = k + 1;
                end
                k = 2;
                for i = (p + 1) * (period_TA - 1) + 1:length(x)
                    derPolynom(index).dAR_pol(k-1) = - x(i);
                    k = k + 1;
                end

            else
                for i = 1:p * (period_TA - 1)
                    derPolynom(index).dT_L.Coefficients{k} = x(i);
                    k = k + 1;
                end
                k = 2;
                for i = p * (period_TA - 1) + 1:length(x)
                    derPolynom(index).dAR_pol(k - 1) = - x(i);
                    k = k + 1;
                end

            end
            if (period_TA == 1)
                derPolynom(index).dT_L = LagOp({1});
            end

        end
        if ne(q,0)
            MA_coeff_init = LagOp([1 MA_coeff]);
            C_L = mtimes(T_L, MA_coeff_init);
            [coefficients_C_L] = toCellArray(C_L);
            coefficients_C_L = cell2mat(coefficients_C_L);

        else

            MA_coeff_init = LagOp([1]);
            [coefficients_C_L] = toCellArray(T_L);
            coefficients_C_L = cell2mat(coefficients_C_L);

        end
        MA_TA_pol = LagOp([1 MA_H]);
        coefficients_MA_TA_pol = cell2mat(toCellArray(MA_TA_pol));

        dBeta = [];
        for i=1:p + q
            x0 = zeros(q_new + 1, 1);

            if i > p
                dMA_L = LagOp({0});

                for j = 1:q+1
                    if j == i - p
                        dMA_L.Coefficients{j} = 1;
                    else
                        dMA_L.Coefficients{j} = 0;
                    end
                end
                dC_L = mtimes(dMA_L, T_L);
                [coefficients_dC_L] = toCellArray(dC_L);
                coefficients_dC_L = cell2mat(coefficients_dC_L) ;
                difer = length(coefficients_C_L) - length(coefficients_dC_L);
                if difer > 0
                    for kk = 1:difer
                        coefficients_dC_L = [coefficients_dC_L 0];
                    end
                else if difer < 0
                        for kk = 1:abs(difer)
                            coefficients_C_L = [coefficients_C_L 0];
                        end
                    end
                end

            else
                dC_L = mtimes(derPolynom(i).dT_L, MA_coeff_init);
                [coefficients_dC_L] = toCellArray(dC_L);
                coefficients_dC_L = cell2mat(coefficients_dC_L);
                difer = length(coefficients_C_L) - length(coefficients_dC_L);
                if difer > 0
                    for kk = 1:difer
                        coefficients_dC_L = [coefficients_dC_L 0];
                    end
                else if difer < 0
                        for kk = 1:abs(difer)
                            coefficients_C_L = [coefficients_C_L 0];
                        end
                    end
                end
            end

            f = @(y)aggregationEquationMA(y, q_new, period_TA, coefficients_C_L, coefficients_dC_L, MA_TA_pol, sigma2, sigma2_H);
            %set up the options for the solver
            options = optimset('Display', 'iter', 'Algorithm', 'trust-region-dogleg', 'TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 200000, 'MaxIter', 200000);
            [y, ~, exitflag, output] = fsolve(f, x0, options);
            
            if exitflag <= 0
                display('problem when solving aggregationEquationMA');
                display(output);
                return;
            end
            dMA_pol = [0];
            k = 1;

            for j = 1:length(y) - 1
                dMA_pol(k) = y(j);
                k = k + 1;
            end
            if i<=p
                temp = -derPolynom(i).dAR_pol(:);
            else
                temp = zeros(1, p);
            end
            dBeta = [dBeta(:); temp(:); dMA_pol(:)];
        end
    end

    

