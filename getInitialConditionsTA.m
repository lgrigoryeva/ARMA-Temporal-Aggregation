function [structInitCond] = getInitialConditionsTA(AR, MA, sigma2, freq_TA, pathLength, aggrType, numTrials)
    
    p = length(AR);
    q = length(MA);
    
    CoeffARMA.AR = AR;
    CoeffARMA.MA = MA;
    CoeffARMA.K = sigma2;
    
    structInitCond = [];
    
    j = 1;
    
    while j <= numTrials

        if strcmp(aggrType, 'flow')
            q_new = idivide(int16((p + 1) * (freq_TA - 1) + q), freq_TA);
        else
            q_new = idivide(int16(p * (freq_TA - 1) + q), freq_TA);
        end
        simSeries = simulateARMA(AR, MA, sigma2, pathLength, 1);
        aggrSimSeries = aggregateData(aggrType, freq_TA, simSeries);
        SpecFit = garchset('VarianceModel', 'Constant', 'R', double(p), 'M', double(q_new), 'Display', 'Off');
        CoeffFit = garchfit(SpecFit, aggrSimSeries(:, 1));
        
        if ~(isfield(CoeffFit, 'MA') && q_new)
            CoeffFit.MA = 0;
        end
        
        if ~(isfield(CoeffFit, 'AR') && p)
            CoeffFit.AR = 0;
        end
        
        [AR_ta, MA_ta] = temporalAggregation(p, q, CoeffARMA, freq_TA, aggrType, CoeffFit);
        
        roots_AR_ta = roots([1; -AR_ta(:)]);
        roots_MA_ta = roots([1; MA_ta(:)]);
        
        if all(abs(roots_AR_ta) < 1) && all(abs(roots_MA_ta) < 1)
            structInitCond.MA = CoeffFit.MA;
            structInitCond.K = CoeffFit.K;
            break;
        else
            if (j == numTrials)
                structInitCond.MA = zeros(1, q_new);
                structInitCond.K = 0;
            else
                pathLength = pathLength * 2;
                structInitCond = [];
            end
            j = j + 1;
        end
        
    end
    