function [Psi_p, Psi_q ] = getAutocovMixedRecur(AR, MA, h, precision, MA_el)
% This function is used in the function getAsymptMLEcovar for computing the
% (1,2),(2,1) blocks of SIGMA_beta in a recursive way, see page 7 of Paper 2

    if MA_el < h
        MA_el = h;
    end

    Psi_p = [1; reshape(garchma(AR, [], MA_el + h - 1), [], 1)];
    Psi_q = [1; reshape(garchma(MA, [], MA_el + h - 1), [], 1)];
    
    autocov_el = Psi_p.*Psi_q;
    
    if abs(autocov_el(end)) > precision
        MA_el = MA_el + 1;
        [Psi_p, Psi_q] = getAutocovMixedRecur(AR, MA, h, precision, MA_el);
    end
    
    return 
    
end