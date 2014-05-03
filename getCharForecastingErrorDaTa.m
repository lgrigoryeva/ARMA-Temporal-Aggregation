function [msfeCharact] = getCharForecastingErrorDaTa(Psi, sigma2, w)
    
    K = length(w);
    
    A_char = zeros(K);
    
    for h = 1:K
        temp_ar = [];
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
    
    msfeCharact = sigma2 * w * A_char * w';
    
end