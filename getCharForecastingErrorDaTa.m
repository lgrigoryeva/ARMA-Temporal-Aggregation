function charactMSFE_AD = getCharForecastingErrorDaTa(Psi, sigma2, w)
% This function evaluates the characteristic error associated to the
% all-disaggregated forecasting scheme (see A^{char}_{hh'} in Equation
% (2.18) of Paper 1)

% w - is the aggregation vector
% sigma2 - unconditional variance of the noise
% Psi - the parameters of the MA-representation of ARMA model
% for stock aggregation w contains zero-valued elements up to index K, 
% which is the aggregation period, and for which w(K)=1
% for flow aggregation w contains K one-valued elements

    % period of aggregation
    K = length(w);
    
    A_char = zeros(K);
    
    for h = 1:K
        temp_ar = zeros(K - h + 1, 1);
        for hh = h:K
            s1 = 0;
            for i = 1:h
                temp = - h + i + hh;
                
                for ii = 1:hh
                    if (ii == temp)
                        s1 = s1 + Psi(ii) * Psi(i);
                    end
                end
                
            end
            temp_ar(hh) = s1;
        end
        A_char(h, :) = temp_ar;
    end
    A_char = (A_char + A_char') - diag(diag(A_char));
    
    charactMSFE_AD = sigma2 * w * A_char * w';
    
end