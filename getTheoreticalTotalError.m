function [totalMSFE_da, totalMSFE_ta, totalMSFE_ta_opt, ...
    charactMSFE_da, charactMSFE_ta, charactMSFE_ta_opt, optSteps] = ...
                getTheoreticalTotalError(AR, MA, sigma2, aggrType, tStart, ...
                                    tEnd, tStep, hStart, hEnd, hStep, pathLength, numTrials, P_cut)    
    startMinValue = 10e6;
    
    p = length(AR);
    q = length(MA);
    r = max(p, q);
    
    
    allowedRootsValue = 0.999;
    checkZeros = 1;
    
    %forecasting horizon
    hNum = (hEnd - hStart + hStep) / hStep;

    %modelled sample length
    tNum = floor((tEnd - tStart + tStep - hEnd) / tStep);

    % precisions and tolerance values
    precConvToMA    = 1.0e-12;
    recursStartLevel= 100;
    
    P_max = tStart + hNum - 1 + r;
    if nargin < 13
        P_cut = tStart + r;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %
    % INITIAL COMPUTATIONS OF THE PARAMETERS, JACOBIANS, ASYMPTOTIC
    % COVARIANCE MATRIX AND COVARIANCE MATRIX FOR DELTA NORMAL
    % METHOD TO BE USED
    %
    
    % canonical asymptotic covariance matrix of MLE estimator
    V_beta = getAsymptMLEcovar(AR, MA, precConvToMA, recursStartLevel);    
    
    [Psi, Pi] = getPsiPi(AR, MA, P_max);
    [psiJac, piJac] = getJacobianPsiPi(AR, MA, Psi, Pi, P_max);

    for i = length(psiJac) + 1:P_max + 1
        psiJac(i, :) = 0;
    end
    for i = length(piJac) + 1:P_max + 1
        piJac(i, :) = 0;
    end
    covMatrixPsiPi = [psiJac; piJac] * V_beta * [psiJac; piJac]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %
    % CHARACTERISTIC FORECASTING AND ESTIMATION ERROR COMPUTATION 
    % FOR THE DISAGGREGATED MODEL AND DATA
    
    %
    % CHARACTERISTIC FORECASTING ERROR
    %
    
    theorMSFE_da = zeros(1, hNum);
    
    for h = 1:hNum
        hLen = hStart + (h - 1) * hStep;
        w = zeros(1, hLen);
        if strcmp(aggrType, 'stock')
            w(1, end) = 1;
        else
            w = ones(1, hLen);
        end
        theorMSFE_da(h) = getCharForecastingErrorDaTa(Psi, sigma2, w);
    end
    
    %
    % ESTIMATION FORECASTING ERROR - GO FORMULA 
    % (LUTKEPOHL WITHOUT MONTE CARLO)
    % 
    
    erCorrec_da_GOlutkepohl = zeros(tNum, hNum);
        
    for t = 1:tNum
        errorEstGOLutkepohl_da = zeros(hNum, 1);
        tLength     = (t - 1) * tStep + tStart;

        for i = 1:hNum
            hLength_ = hStart + (i - 1) * hStep;
            w = zeros(1, hLength_);
            if strcmp(aggrType, 'stock')
                w(1, end) = 1;
            else
                w = ones(1, hLength_);
            end
            %computation of th etotal MSFE error due to Grigoryeva&Ortega paper
            errorEstGOLutkepohl_da(i) = getEstimationErrorDaTa(Psi, sigma2, covMatrixPsiPi, w, P_cut, tLength, P_max);
        end
        disp(t);
        erCorrec_da_GOlutkepohl(t, :) = errorEstGOLutkepohl_da';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %
    % MODEL AGGREGATION BLOCK FOR GDP DATA
    %

    % disaggregated model parameters
    paramDA.AR = AR;
    paramDA.MA = MA;
    paramDA.K = sigma2;

    dBetaTotal  = struct('dBeta', []);
    paramTA = struct('AR', [], 'MA', [], 'p', [], 'q', [], 'Psi', [], 'Pi', [], 'sigma2', []);

    for h = 1:hNum
        dBetaTotal(h) = struct('dBeta', []);
        paramTA(h) = struct('AR', [], 'MA', [], 'p', [], 'q', [], 'Psi', [], 'Pi', [], 'sigma2', []);
        
        hLen = hStart + (h - 1) * hStep;
        structInitCond = getInitialConditionsTA(paramDA.AR, paramDA.MA, paramDA.K, hLen, pathLength, aggrType, numTrials);

        
        [paramTA(h).AR, paramTA(h).MA, paramTA(h).sigma2, ~, dBeta] = temporalAggregation(p, q, paramDA, hLen, aggrType, structInitCond);
        if paramTA(h).AR == 0
            
            paramTA(h).AR = [];
        end
        
        if checkCoeffRoots(paramTA(h).AR, paramTA(h).MA, allowedRootsValue, checkZeros)
            if isempty(paramTA(h).MA(:))
                paramTA(h).MA(:) = 0;
            end
        end
        [paramTA(h).Psi, paramTA(h).Pi] = getPsiPi(paramTA(h).AR, paramTA(h).MA, P_max);
        dBetaTotal(h).dBeta = dBeta(:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %
    % ESTIMATION ERROR COMPUTATION FOR THE AGGREGATED MODEL AND DATA
    %
    
    errorEstGOLutkepohl_ta = zeros(hNum, 1);
    erCorrec_ta_GOlutkepohl = zeros(tNum, hNum);

    for t = 1:tNum
        tLength = (t - 1) * tStep + tStart;
        for h = 1:hNum

            MA_ = paramTA(h).MA(:);

            if ~p || (p == 1 && paramTA(h).AR == 0)
                AR_ = [];
            else
                AR_ = paramTA(h).AR(:);
            end

            p_est = length(AR_);
            q_est = length(MA_);
            dBeta_reshaped = reshape(dBetaTotal(h).dBeta, p_est + q_est, p + q);
            V_beta_ta = dBeta_reshaped * V_beta * dBeta_reshaped';
    
            [psiJac, piJac] = getJacobianPsiPi(AR_, MA_, paramTA(h).Psi, paramTA(h).Pi, P_max);

            for i = length(psiJac) + 1:P_max + 1
                psiJac(i, :) = 0;
            end
            for i = length(piJac) + 1:P_max + 1
                piJac(i, :) = 0;
            end
            covMatrixPsiPi_ = [psiJac; piJac] * V_beta_ta * [psiJac; piJac]';
             % computation of the total MSFE error due to Grigoryeva&Ortega paper
            errorEstGOLutkepohl_ta(h) = getEstimationErrorDaTa(paramTA(h).Psi, paramTA(h).sigma2, covMatrixPsiPi_, 1, P_cut, tLength, P_max);

            disp(h);
        end

        erCorrec_ta_GOlutkepohl(t, :) = errorEstGOLutkepohl_ta';
    end
    
    disp('estimation error computation is done');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %
    % THEORETICAL MSFE COMPUTATION
    % FORECASTING OF AGGREGATES FOR AGGREGATED MODEL
    % FINITE SAMPLE
    %

    theorMSFE_ta = zeros(hNum, 1);

    for h = 1:hNum
        if checkCoeffRoots(paramTA(h).AR, paramTA(h).MA, allowedRootsValue, checkZeros)
            theorMSFE_ta(h) = getARMAmsfeFinSample(1, paramTA(h).AR(:), paramTA(h).MA(:), paramTA(h).sigma2);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %
    % THEORETICAL MSFE COMPUTATION
    % FORECASTING OF AGGREGATES FOR AGGREGATED MODEL
    % FINITE SAMPLE
    % SELECTION OF THE BEST AGGREGATION SCHEME IN TERMS OF THE STEP/AGGREGATION
    % PERIOD VALUE IF H STEP HAS A FEW MULTIPLIERS, CALCULATION OF THE MINIMAL
    % ERROR VALUE AND OPTIMAL AGGREGATION PERIOD LENGTH
    %
    
    multArr = [2, 3, 4, 5, 6, 7, 8, 9, 10];
    multArrLength = length(multArr);
    erCorrec_ta_opt_GOlutkepohl = zeros(tNum, 1);
    msfError = zeros(hNum, 1);
    curParamTA = struct('AR', [], 'MA', [], 'p', [], 'q', [], 'Psi', [], 'Pi', [], 'sigma2', []);
    MSFE_ta_opt = struct('msfError', msfError, 'erCorrec_ta_opt_GOlutkepohl', erCorrec_ta_opt_GOlutkepohl, 'erCorrec_ta_opt_GOlutkepohl_temp', []);
    tempAggrStruct = struct('stepValues', [], 'MSFE_ta_opt', MSFE_ta_opt, 'curParamTA', curParamTA);

    for h = 1:hNum
        msfError = zeros(hNum, 1);
        curParamTA = struct('AR', [], 'MA', [], 'p', [], 'q', [], 'Psi', [], 'Pi', [], 'sigma2', []);
        MSFE_ta_opt = struct('msfError', msfError, 'erCorrec_ta_opt_GOlutkepohl', erCorrec_ta_opt_GOlutkepohl, 'erCorrec_ta_opt_GOlutkepohl_temp', []);
        tempAggrStruct(h) = struct('stepValues', [], 'MSFE_ta_opt', MSFE_ta_opt, 'curParamTA', curParamTA);
        
        hLengthCur = hStart + (h - 1) * hStep;
        tempAggrStruct(h).MSFE_ta_opt = [];
        modArr = mod(hLengthCur, multArr);

        if (hLengthCur == 1)
            tempAggrStruct(h).stepValues = hLengthCur;
        else
            for i = 1:multArrLength
                if ~modArr(i)
                    tempAggrStruct(h).stepValues = [tempAggrStruct(h).stepValues, multArr(i)];
                end
            end
        end
        
        numVariants = length(tempAggrStruct(h).stepValues);
        
        for i = 1:numVariants
            curStep = tempAggrStruct(h).stepValues(i);
            h_div = hLengthCur/curStep;
            w = zeros(1, h_div);
            if strcmp(aggrType, 'stock')
                w(1, end) = 1;
            else
                w = ones(1, h_div);
            end
            tempAggrStruct(h).MSFE_ta_opt(curStep).msfError = getCharForecastingErrorDaTa(paramTA(curStep).Psi, paramTA(curStep).sigma2, w);
            tempAggrStruct(h).curParamTA(curStep) = paramTA(curStep);
            
            if (length(tempAggrStruct(h).curParamTA(curStep).AR) == 1 && tempAggrStruct(h).curParamTA(curStep).AR == 0)
                tempAggrStruct(h).curParamTA(curStep).p = 0;
                tempAggrStruct(h).curParamTA(curStep).AR = [];

            else
                tempAggrStruct(h).curParamTA(curStep).p = length(tempAggrStruct(h).curParamTA(curStep).AR);
            end
            
            tempAggrStruct(h).curParamTA(curStep).q = length(tempAggrStruct(h).curParamTA(curStep).MA);
        end
        
    end

    for i = 1:hNum
        hLength = hStart + (i - 1) * hStep;
        curTempAggrStruct = tempAggrStruct(i);
        
        
        for j = 1:length(curTempAggrStruct.stepValues)
            curStep = curTempAggrStruct.stepValues(j);
            hLengthCur = int8(hLength/curStep);
            curParamTA_ = curTempAggrStruct.curParamTA(curStep);
            w = zeros(1, hLengthCur);
            if strcmp(aggrType, 'stock')
                w(1, end) = 1;
            else
                w = ones(1, hLengthCur);
            end
            curdBeta_reshaped = reshape(dBetaTotal(curStep).dBeta, curParamTA_.p + curParamTA_.q, p + q);
            curV_beta_ta = curdBeta_reshaped * V_beta * curdBeta_reshaped';
            
            [curPsiJac, curPiJac] = getJacobianPsiPi(curParamTA_.AR, curParamTA_.MA, curParamTA_.Psi, curParamTA_.Pi, P_max);
            
            for m = length(curPsiJac) + 1:P_max + 1
                curPsiJac(m, :) = 0;
            end
            for m = length(curPiJac) + 1:P_max + 1
                curPiJac(m, :) = 0;
            end
            curCovMatrixPsiPi_ = [curPsiJac; curPiJac] * curV_beta_ta * [curPsiJac; curPiJac]';
            
            % computation of the total MSFE error due to Grigoryeva&Ortega paper
            tempAggrStruct(i).MSFE_ta_opt(curStep).erCorrec_ta_opt_GOlutkepohl_temp = getEstimationErrorDaTa(curParamTA_.Psi, curParamTA_.sigma2, curCovMatrixPsiPi_, w, P_cut, 1, P_max);
        end
    end
    
    totalOptError = zeros(tNum, hNum);
    optSteps = zeros(hNum, 1);
    theorMSFE_ta_opt = zeros(tNum, hNum);
    
    for t = 1:tNum
        tLength = (t - 1) * tStep + tStart;
        for h = 1:hNum
            curTempAggrStruct = tempAggrStruct(h);
            totalOptError(t, h) = startMinValue;
            
            for j = 1:length(curTempAggrStruct.stepValues)
                curStep = curTempAggrStruct.stepValues(j);
                if (totalOptError(t, h) > (tempAggrStruct(h).MSFE_ta_opt(curStep).erCorrec_ta_opt_GOlutkepohl_temp / tLength + tempAggrStruct(h).MSFE_ta_opt(curStep).msfError))
                   totalOptError(t, h) = tempAggrStruct(h).MSFE_ta_opt(curStep).erCorrec_ta_opt_GOlutkepohl_temp / tLength  + tempAggrStruct(h).MSFE_ta_opt(curStep).msfError;
                   optSteps(h) = curStep;
                end
            end
            
            theorMSFE_ta_opt(t, h) = tempAggrStruct(h).MSFE_ta_opt(optSteps(h)).msfError;
        end
    end
    
    totalMSFE_ta = repmat(theorMSFE_ta', tNum, 1) + erCorrec_ta_GOlutkepohl;
    totalMSFE_da = repmat(theorMSFE_da, tNum, 1) + erCorrec_da_GOlutkepohl;
    totalMSFE_ta_opt = totalOptError;
    charactMSFE_ta = theorMSFE_ta';
    charactMSFE_da = theorMSFE_da;
    charactMSFE_ta_opt = theorMSFE_ta_opt;
end

    
