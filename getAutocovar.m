function  autocov  = getAutocovar(AR, MA, h, precision, recursStartLevel)
% Provides the autocovariance function of an ARMA model up to lag h
% precision sets the tolerance accepted in the use of the AR infinity representation
    if recursStartLevel < h
        recursStartLevel = h;
    end
    Psi = [1; reshape(garchma(AR, MA, recursStartLevel + h - 1), [], 1)];

    autocov_el = Psi(1:recursStartLevel).*Psi(1 + h:end);
    
    if abs(autocov_el(end)) > precision
        recursStartLevel = recursStartLevel + 1;
        autocov = getAutocovar(AR, MA, h, precision, recursStartLevel);
    else
        autocov = sum(autocov_el);
    end
    return 
end

