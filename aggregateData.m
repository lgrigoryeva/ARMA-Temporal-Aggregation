function [aggr_sample] = aggregateData(aggr_type, period_TA, sample)
% stock and flow aggregation of a sample with a prescribed period

    if (period_TA == 1)
        aggr_sample = sample;
        return;
    end
        
        
    [T, N] = size(sample);
    sample = flipud(sample);
    
    if strcmp('stock', aggr_type)
        aggr_sample = zeros(idivide(int16(T - 1), period_TA) + 1, N);
        for ii = 1:N
            
            
            aggr_sample(1, ii) = sample(1, ii);
            j = 2;
            for i = period_TA:T
                
                if ~mod(i - 1, period_TA)
                    aggr_sample(j, ii) = sample(i, ii);
                    j = j + 1;
                end
            end
        end
        aggr_sample = flipud(aggr_sample);
    else
        aggr_sample = zeros(idivide(int16(T), period_TA), N);
        for ii = 1:N
            aggr = 0;
            k = 0;
            j = 1;
            for i = 1:T
                aggr = aggr + sample(i, ii);
                k = k + 1;
                if (k == period_TA)
                    aggr_sample(j, ii) = aggr;
                    j = j + 1;
                    aggr = 0;
                    k = 0;
                end
            end
       end
        aggr_sample = flipud(aggr_sample);
    end
end

