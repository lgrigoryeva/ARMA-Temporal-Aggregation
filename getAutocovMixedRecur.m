function [Psi_p, Psi_q ] = getAutocovMixedRecur(AR, MA, h, precision, recursStartLevel)
% This function is used in the function getAsymptMLEcovar for computing the
% (1,2),(2,1) blocks of SIGMA_beta in a recursive way, see page 7 of Paper 2

    if recursStartLevel < h
        recursStartLevel = h;
    end

    Psi_p = [1; reshape(garchma(AR, [], recursStartLevel + h - 1), [], 1)];
    Psi_q = [1; reshape(garchma(MA, [], recursStartLevel + h - 1), [], 1)];
    
    autocov_el = Psi_p.*Psi_q;
    
    if abs(autocov_el(end)) > precision
        recursStartLevel = recursStartLevel + 1;
        [Psi_p, Psi_q] = getAutocovMixedRecur(AR, MA, h, precision, recursStartLevel);
    end
    
    return 
    
end