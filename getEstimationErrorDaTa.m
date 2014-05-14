function estimMSFE = getEstimationErrorDaTa(Psi, sigma2, covPsiPiMatrix, w, P_cut, T, P_max)
% This function evaluates the estimation error associated to the
% (1)all-disaggregated, (2)all-aggregated, or (3)hybrid forecasting schemes, respectively, 
% for (1): see Equation (2.18) of Paper 1 for D, F, G in the body of this function;
% for (2): see square bracket in Equation (2.20) of Paper 1;
% for (3): see square bracket in Equation (2.13) of Paper 1.

% w - is the aggregation vector
% for (1):
% in the case of stock aggregation w contains zero-valued elements up to index K, 
% which is the aggregation period, and for which w(K) = 1
% in the case of flow aggregation w contains K one-valued elements;
% for (2) and (3): 
% set w = 1.

% sigma2 - unconditional variance of the noise
% Psi - the parameters of the MA-representation of considered ARMA model
    Psi = Psi(2:end);
    
    K = length(w);
    D = zeros(K);
    for h = 1:K
        P = P_cut + h - 1;
        temp_ar = zeros(1, K);
        for hh = h:K
            PP = P_cut + hh - 1;
            
            s1 = 0;
            for i = h:P
                
                for ii = hh:PP
                    if (h - i == hh - ii)
                        s1 = s1 + covPsiPiMatrix(i + 1, ii + 1);
                    end
                end
            end
            temp_ar(hh) = s1;
        end
        D(h, :) = temp_ar;
    end
    D = 1/T * ((D + D') - diag(diag(D)));
    
    
    F = zeros(K);    
    for h = 1:K
        P = P_cut + h - 1;
        for hh = 1:K
            PP = P_cut + hh - 1;
            s1 = 0;
            for i = h:P
                covCols = covPsiPiMatrix(i + 1, :);
                for ii = hh:PP
                    for jj = 0:PP - ii
                        for kk = 0:PP - ii - jj
                            
                            if (h - i == hh - ii - jj - kk)
                                if ~kk
                                    s1 = s1 + covCols(:, P_max + jj + 2) * Psi(ii);
                                else
                                    s1 = s1 + Psi(kk) * covCols(:, P_max + jj + 2) * Psi(ii);
                                end
                            end
                        end
                    end
                end
            end
            F(h, hh) = s1;
        end
    end    
    F = F * 2 / T;
    
    G = zeros(K);
    for h = 1:K
        P = P_cut + h - 1;
        temp_ar = zeros(1, K);

        for hh = h:K
            PP = P_cut + hh - 1;
            s4 = 0;
            for i = h:P
                for j = 0:P - i
                    for k = 0:P - i - j
                        for ii = hh:PP
                            for jj = 0:PP - ii
                                for kk = 0:PP - ii - jj
                                    if (h - i - j - k == hh - ii - jj - kk)
                                        if ~k
                                            Psi_k = 1;
                                        else
                                            Psi_k = Psi(k);
                                        end
                                        if ~kk
                                            Psi_kk = 1;
                                        else
                                            Psi_kk = Psi(kk);
                                        end
                                        s4 = s4 + Psi(i) * Psi_k * Psi(ii) * Psi_kk * covPsiPiMatrix(j + P_max + 2, jj + P_max + 2);
                                    end
                                end
                            end
                        end
                    end
                end
            end
            temp_ar(hh) = s4;
        end
        G(h, :) = temp_ar;
    end
    G = 1 / T * ((G + G') - diag(diag(G)));

    estimMSFE = sigma2 * w * (F + D + G) * w';
        
end