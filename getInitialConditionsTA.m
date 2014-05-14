function structInitCond = getInitialConditionsTA(AR, MA, sigma2, period_TA, pathLength, aggrType, numTrials)
    
% This function generates initial values for the nonlinear system which
% has to be solved in order to obtain the coefficients of the new
% temporally aggregated ARMA model (see Equation (2.14) of Paper 2 for example)

% In order to obtain the initial values we generate a sample of length
% pathLength, aggregate it according to a given aggregation type "stock" or
% "flow" given in aggr_type, and subsequently fit an ARMA(p,q*) model to
% it, where q* is the order of MA part of the new aggregated model.

% If the outcome of this procedure gives the ARMA(p,q*) model which is not
% causal or invertible, we repeat this procedure by enlarging pathLength by
% pathLength=pathLength*2 and as many times as given by numTrials

    p = length(AR);
    q = length(MA);
    
    CoeffARMA.AR = AR;
    CoeffARMA.MA = MA;
    CoeffARMA.K = sigma2;
    
    structInitCond = [];
    
    j = 1;
    
    while j <= numTrials

        if strcmp(aggrType, 'flow')
            q_new = idivide(int16((p + 1) * (period_TA - 1) + q), period_TA);
        else
            q_new = idivide(int16(p * (period_TA - 1) + q), period_TA);
        end
        simSeries = simulateARMA(AR, MA, sigma2, pathLength, 1);
        aggrSimSeries = aggregateData(aggrType, period_TA, simSeries);
        SpecFit = garchset('VarianceModel', 'Constant', 'R', double(p), 'M', double(q_new), 'Display', 'Off');
        CoeffFit = garchfit(SpecFit, aggrSimSeries(:, 1));
        
        if ~(isfield(CoeffFit, 'MA') && q_new)
            CoeffFit.MA = 0;
        end
        
        if ~(isfield(CoeffFit, 'AR') && p)
            CoeffFit.AR = 0;
        end
        if strcmp(aggrType, 'stock')
            w = zeros(period_TA, 1);
            w(end) = 1;
        else
            w = ones(period_TA, 1);
        end
        [AR_ta, MA_ta] = temporalAggregation(p, q, CoeffARMA, w, CoeffFit);
        
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
    