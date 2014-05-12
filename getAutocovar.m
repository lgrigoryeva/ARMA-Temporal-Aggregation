function  autocov  = getAutocovar(AR, MA, h, precision, MA_el)
% Provides the autocovariance function of an ARMA model up to lag h
% precision sets the tolerance accepted in the use of the AR infinity representation
    if MA_el < h
        MA_el = h;
    end
    Psi = [1; reshape(garchma(AR, MA, MA_el + h - 1), [], 1)];

    autocov_el = Psi(1:MA_el).*Psi(1 + h:end);
    
    if abs(autocov_el(end)) > precision
        MA_el = MA_el + 1;
        autocov = getAutocovar(AR, MA, h, precision, MA_el);
    else
        autocov = sum(autocov_el);
    end
    return 
end

