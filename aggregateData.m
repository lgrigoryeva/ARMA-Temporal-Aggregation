function [X_initial_aggr] = aggregateData(sampl_type, freq_TA, X_init)
    if (freq_TA == 1)
        X_initial_aggr = X_init;
        return;
    end
        
        
    [T, N] = size(X_init);
    X_init = flipud(X_init);
    
    if strcmp('stock', sampl_type)
        X_initial_aggr = zeros(idivide(int16(T - 1), freq_TA) + 1, N);
        for ii = 1:N
            
            
            X_initial_aggr(1, ii) = X_init(1, ii);
            j = 2;
            for i = freq_TA:T
                
                if ~mod(i - 1, freq_TA)
                    X_initial_aggr(j, ii) = X_init(i, ii);
                    j = j + 1;
                end
            end
        end
        X_initial_aggr = flipud(X_initial_aggr);
    else
        X_initial_aggr = zeros(idivide(int16(T), freq_TA), N);
        for ii = 1:N
            aggr = 0;
            k = 0;
            j = 1;
            for i = 1:T
                aggr = aggr + X_init(i, ii);
                k = k + 1;
                if (k == freq_TA)
                    X_initial_aggr(j, ii) = aggr;
                    j = j + 1;
                    aggr = 0;
                    k = 0;
                end
            end
       end
        X_initial_aggr = flipud(X_initial_aggr);
    end
end

