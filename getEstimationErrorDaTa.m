function [msfeEstim] = getEstimationErrorDaTa(Psi, sigma2, covPsiPiMatrix, w, P_cut, T_sample, P_max)
    

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
    D = 1/T_sample * ((D + D') - diag(diag(D)));
    
    
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
    F = F * 2 / T_sample;
    
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
    G = 1 / T_sample * ((G + G') - diag(diag(G)));

    msfeEstim = sigma2 * w * (F + D + G) * w';
        
end