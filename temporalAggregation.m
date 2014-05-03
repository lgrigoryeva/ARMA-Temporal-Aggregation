function [AR_pol, MA_pol, var_e_TA, T_L, dBeta] = temporalAggregation(p, q, Coeff, freq_TA, aggrType, structInitCond)
    %
    %PARAMETERS INITIALIZATION FOR TEMPORARY AGGREGATION PROCEDURE
    %
    %take the parameters from the fitted ones
    AR_coeff = -Coeff.AR;
    if isfield(Coeff, 'MA') && q
        MA_coeff = Coeff.MA;
    else
        MA_coeff=0;
    end
    var_e = Coeff.K;

    %
    %TEMPORAL AGGREGATION PROCEDURE
    %
    %set initial point
    if strcmp(aggrType, 'flow')
         if (((p + 1) * (freq_TA - 1)) + p == 0)
            x0 = 0;
         else
            x0 = -1 * ones(((p + 1) * (freq_TA - 1)) + p, 1);
         end
    else
        if ~p
            x0 = 0;
        else
            x0 = -10 * ones((p * (freq_TA - 1)) + p, 1);
        end
    end
    
    %get the T(L) aggregation lag operator and A*(L)=T(L)*A(L) lag operator for aggregated AR
    %create auxiliary function
    f = @(x)getAggrPolynoms(x, AR_coeff, freq_TA, aggrType);
    %set up options for the solver
    options = optimset('Display', 'off', 'Algorithm', 'levenberg-marquardt', 'TolFun', 1e-12, 'TolX', 1e-12);
    %get the solution for the coefficients
    [x, ~, exitflag, output] = fsolve(f, x0, options);
    if ~exitflag
        display('problem when solving getAggrPolynoms');
        display(output);
        return;
    end

    %
    %PROCESSING THE RESULT x
    %
    T_L = LagOp({1});
    coef_TA(1) = 1;
    AR_pol = [0];
    k = 1;
    if strcmp(aggrType, 'flow')
        for i = 1:(p + 1) * (freq_TA - 1)
            T_L.Coefficients{k} = x(i);
            coef_TA(k + 1) = x(i);
            k = k + 1;
        end
        k = 2;
        for i = (p + 1) * (freq_TA - 1) + 1:length(x)
            AR_pol(k-1) = - x(i);
            k = k + 1;
        end
        %new order of temporarily aggregated ARMA
        p_new = p;
        q_new = idivide(int16((p + 1) * (freq_TA - 1) + q), freq_TA);
    else
        for i = 1:p * (freq_TA - 1)
            T_L.Coefficients{k} = x(i);
            coef_TA(k + 1) = x(i);
            k = k + 1;
        end
        k = 2;
        for i = p * (freq_TA - 1) + 1:length(x)
            AR_pol(k - 1) = - x(i);
            k = k + 1;
        end
        %new order of temporarily aggregated ARMA
        p_new = p;
        q_new = idivide(int16(p * (freq_TA - 1) + q), freq_TA);
    end


    %get the new MA aggregated polynomial
    if ne(q,0)
        MA_L = mtimes(T_L, LagOp([1 MA_coeff]));
        [coefficients_TA] = toCellArray(MA_L);
    else
        [coefficients_TA] = toCellArray(T_L);
    end

    %
    %TEMPORAL AGGREGATION LOWERING THE ORDER OF L^k = B
    %
    %define the variance of the TA model
    coefficients_TA = cell2mat(coefficients_TA);
    %set up initial point for the solver
    
    x0 = [structInitCond.MA, structInitCond.K];
    %set the auxiliary function for estimation of the parameters of the new
    %model
    f = @(y)aggregationEquation(y, q_new, coefficients_TA, var_e, freq_TA);
    %set up the options for the solver
    options = optimset('Display', 'off', 'Algorithm', 'levenberg-marquardt', 'TolFun', 1e-20, 'TolX', 1e-20);
    [y, ~, exitflag, output] = fsolve(f, x0, options);
    if ~exitflag
        display('problem when solving aggregationEquation');
        display(output);
        return;
    end
    %y contains parameters of the aggregated model and variance of the
    %aggregated residuals

    MA_pol = [0];
    k = 1;
    for i = 1:q_new
        MA_pol(k) = (y(i));
        k = k + 1;
    end

    var_e_TA = y(q_new + 1);

    if nargout > 4
        if strcmp(aggrType, 'flow')
            x0 = -10 * ones(((p + 1) * (freq_TA - 1)) + p, 1);
        else
            x0 = -10 * ones((p * (freq_TA - 1)) + p, 1);
        end
        derPolynom = [];
        for index = 1:p

            f = @(x)getAggrPolynomsDeriv(x, AR_coeff, coef_TA, freq_TA, aggrType, index);
            options = optimset('Display', 'off', 'Algorithm', 'levenberg-marquardt', 'TolFun', 1e-12, 'TolX', 1e-12);
            [x, ~, exitflag, output] = fsolve(f, x0, options);
            
            if ~exitflag
                display('problem when solving getAggrPolynomsDeriv');
                display(output);
                return;
            end
            
            
            derPolynom(index).dT_L = LagOp({0});
            derPolynom(index).dAR_pol = [0];
            k = 1;
            if strcmp(aggrType, 'flow')
                for i = 1:(p + 1) * (freq_TA - 1)
                    derPolynom(index).dT_L.Coefficients{k} = x(i);
                    k = k + 1;
                end
                k = 2;
                for i = (p + 1) * (freq_TA - 1) + 1:length(x)
                    derPolynom(index).dAR_pol(k-1) = - x(i);
                    k = k + 1;
                end

            else
                for i = 1:p * (freq_TA - 1)
                    derPolynom(index).dT_L.Coefficients{k} = x(i);
                    k = k + 1;
                end
                k = 2;
                for i = p * (freq_TA - 1) + 1:length(x)
                    derPolynom(index).dAR_pol(k - 1) = - x(i);
                    k = k + 1;
                end
                %new order of temporarily aggregated ARMA

            end
            if (freq_TA == 1)
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
        MA_TA_pol = LagOp([1 MA_pol]);
        coefficients_MA_TA_pol = cell2mat(toCellArray(MA_TA_pol));

        %dBeta = [-dAR_pol(:)];
        dBeta = [];
        for i=1:p + q
            x0 = zeros(q_new + 1, 1);
            %set the auxiliary function for estimation of the parameters of the new
            %model

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

            f = @(y)aggregationEquationMA(y, q_new, freq_TA, coefficients_C_L, coefficients_dC_L, MA_TA_pol, var_e, var_e_TA);
            %set up the options for the solver
            options = optimset('Display', 'off', 'LargeScale','off', 'Algorithm', 'levenberg-marquardt');
            [y, ~, exitflag, output] = fsolve(f, x0, options);
            
            if ~exitflag
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

    

