function [totalMSFE_AD, totalMSFE_H, totalMSFE_OH, ...
    charactMSFE_AD, charactMSFE_H, charactMSFE_OH, optSteps, params_H] = ...
                getTheoreticalTotalError(AR, MA, sigma2, aggrType, tStart, ...
                                    tEnd, tStep, hStart, hEnd, hStep, pathLength, numTrials, P_cut)    
    
% This function evaluats the total and characteristic errors associated to the following forecasting schemes:
% (1) All-disaggregated (AD)
% (2) Hybrid (H)
% (3) Optimal Hybrid (OH)
% Output values: 
% totalMSFE_AD - total forecasting error associated to (1)
% totalMSFE_H - total forecasting error associated to (2)
% totalMSFE_OH - total forecasting error associated to (3)
% charactMSFE_AD - characteristic forecasting error associated to (1)
% charactMSFE_H - characteristic forecasting error associated to (2)
% charactMSFE_OH - characteristic forecasting error associated to (3)
% REMARK: estimation errors associated to (1), (2) and (3), respectively
% can be obtained as the difference between the corresponding total and
% characteristic errors


% Input parameters:
% P_cut - if given, sets the cut-off for the number of Psi and Pi elements
% that are considered at the time of evaluation (this is helpful when the estimation 
% sample size is above 70)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%                    INITIALIZATION BLOCK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
    % determine the order of ARMA
    p = length(AR);
    q = length(MA);
    r = max(p, q);
    
    % upper bound of allowed values for the roots of AR and MA polynomials
    % that guarantee the causality and invertibility of the ARMA model
    allowedRootsValue = 0.999;
    checkZeros = 1;
    
    % length of the set of forecasting horizons of interest
    hNum = (hEnd - hStart + hStep) / hStep;

    % length of the set of estimation sample lengths of interest
    tNum = floor((tEnd - tStart + tStep - hEnd) / tStep);

    % utility precisions and tolerance values (for example for the function getAutocovMixedRecur)
    precConvToMA = 1.0e-12;
    recursStartLevel = 100;
    
    % P_max sets the maximum number of Psi and Pi components to be used
    % throughout the computations (for maxTlength as maximum estimation sample size 
    % of interest and hNum - maximum forecasting horizon of interest)
    maxTlength = (tNum - 1) * tStep + tStart;
    P_max = maxTlength + hNum - 1 + r;
    
    % if P_cut is not given as an input parameter, it is set equal to P_max
    if nargin < 13
        P_cut = P_max;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%                    COMMON INITIAL COMPUTATION BLOCK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % INITIAL COMPUTATIONS OF THE PSI AND PI PARAMETERS, ASSOCIATED JACOBIANS, 
    % ASYMPTOTIC COVARIANCE MATRIX AND COVARIANCE MATRIX FOR DELTA NORMAL
    % METHOD TO BE USED (ALL FOR INITIAL DISAGGREGATED ARMA MODEL)
    
    % canonical asymptotic covariance matrix of MLE estimator
    V_beta = getAsymptMLEcovar(AR, MA, precConvToMA, recursStartLevel);    
    
    % Psi and Pi parameters
    [Psi, Pi] = getPsiPi(AR, MA, P_max);
    
    % Jacobian for Psi and Pi parameters (see J_{XI_P} in Definition (4.8) of Paper 1)
    [psiJac, piJac] = getJacobianPsiPi(AR, MA, Psi, Pi, P_max);

    % if the lengths of psiJac and piJac are smaller than P_max,
    % set all the rest equal to 0
    for i = length(psiJac) + 1:P_max + 1
        psiJac(i, :) = 0;
    end
    for i = length(piJac) + 1:P_max + 1
        piJac(i, :) = 0;
    end
    
    % covariance matrix (see Definition (4.8) of Paper 1)
    covMatrixPsiPi = [psiJac; piJac] * V_beta * [psiJac; piJac]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%                    COMPUTATION OF CHARACTERISTIC AND ESTIMATION ERRORS
    %%%%                    ASSOCIATED TO THE ALL-DISAGGREGATED FORECASTING SCHEME
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    % CHARACTERISTIC FORECASTING ERROR
    % see A^{char}_{hh'} in Equation (2.18) of Paper 1)
    
    charactMSFE_AD = zeros(1, hNum);
    
    for h = 1:hNum
        hLen = hStart + (h - 1) * hStep;
        w = zeros(1, hLen);
        if strcmp(aggrType, 'stock')
            w(1, end) = 1;
        else
            % if aggregation scheme is 'flow'
            w = ones(1, hLen);
        end
        charactMSFE_AD(h) = getCharForecastingErrorDaTa(Psi, sigma2, w);
    end
    
    %
    % ESTIMATION ERROR ASSOCIATED TO THE ALL-DISAGGREGATED FORECASTING SCHEME
    % (Equation (2.18) of Paper 1)
    % 
    
    estimMSFE_AD = zeros(tNum, hNum);
        
    parfor t = 1:tNum
        estimMSFE_AD_h = zeros(hNum, 1);
        tLength = (t - 1) * tStep + tStart;

        for i = 1:hNum
            hLength_ = hStart + (i - 1) * hStep;
            w = zeros(1, hLength_);
            if strcmp(aggrType, 'stock')
                w(1, end) = 1;
            else
                w = ones(1, hLength_);
            end
            estimMSFE_AD_h(i) = getEstimationErrorDaTa(Psi, sigma2, covMatrixPsiPi, w, P_cut, tLength, P_max);
        end
        disp(t);
        estimMSFE_AD(t, :) = estimMSFE_AD_h';
    end
    
    disp('estimation error computation for the all-disaggregated scheme is completed');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%                    BLOCK FOR TEMPORAL AGGREGATION OF THE INITIAL ARMA(p,q) MODEL 
    %%%%                    FOR EACH FORECASTING HORIZON h WE DETERMINE THE
    %%%%                    AGGREGATED MODEL ARMA(p,q*) WITH THE
    %%%%                    CORRESPONDING AR, MA, Psi AND Pi PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initial ARMA model parameters
    param_AD.AR = AR;
    param_AD.MA = MA;
    param_AD.K = sigma2;

    dBetaTotal  = struct('dBeta', []);
    params_H = struct('AR', [], 'MA', [], 'p', [], 'q', [], 'Psi', [], 'Pi', [], 'sigma2', []);

    for h = 1:hNum
        dBetaTotal(h) = struct('dBeta', []);
        params_H(h) = struct('AR', [], 'MA', [], 'p', [], 'q', [], 'Psi', [], 'Pi', [], 'sigma2', []);
        
        hLen = hStart + (h - 1) * hStep;
        if strcmp(aggrType, 'stock')
            w = zeros(hLen, 1);
            w(end) = 1;
        else
            w = ones(hLen, 1);
        end
        % get initial values for the nonlinear temporal aggregation
        % procedure
        structInitCond = getInitialConditionsTA(param_AD.AR, param_AD.MA, param_AD.K, hLen, pathLength, aggrType, numTrials);

        % temporal aggregation procedure
        % dBeta contains Jacobian of temporally aggregated parameters (see J_{beta_Y} in point (ii) of statement of Theorem 2.2 of Paper 1)
        
        [params_H(h).AR, params_H(h).MA, params_H(h).sigma2, ~, dBeta] = temporalAggregation(p, q, param_AD, w, structInitCond);
        if params_H(h).AR == 0
            
            params_H(h).AR = [];
        end
        
        if checkCoeffRoots(params_H(h).AR, params_H(h).MA, allowedRootsValue, checkZeros)
            if isempty(params_H(h).MA(:))
                params_H(h).MA(:) = 0;
            end
        end
        % P_max_h is defined as P in point (ii)  of statement of Theorem 2.2 of Paper 1)
        r_star = max([length(params_H(h).AR), length(params_H(h).MA)]);
        P_max_h = floor(((tNum - 1) * tStep + tStart)/h) + r_star;
        [params_H(h).Psi, params_H(h).Pi] = getPsiPi(params_H(h).AR, params_H(h).MA, P_max_h);
        dBetaTotal(h).dBeta = dBeta(:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%                    COMPUTATION OF CHARACTERISTIC AND ESTIMATION ERRORS
    %%%%                    ASSOCIATED TO THE HYBRID FORECASTING SCHEME
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % CHARACTERISTIC FORECASTING ERROR
    % (see sigma*^2 in Equation (2.13) of Paper 1)
    
    charactMSFE_H = zeros(hNum, 1);

    for h = 1:hNum
        if checkCoeffRoots(params_H(h).AR, params_H(h).MA, allowedRootsValue, checkZeros)
            charactMSFE_H(h) = params_H(h).sigma2;
        else
            display('Problems with the temporal aggregation procedure are encountered');
            totalMSFE_H = [];
            totalMSFE_OH = [];
            charactMSFE_H = [];
            charactMSFE_OH = [];
            optSteps = [];
            return;
        end
    end

    %
    % ESTIMATION ERROR ASSOCIATED TO THE HYBRID FORECASTING SCHEME
    % (see square bracket in Equation (2.13) of Paper 1)
    % 
    
    
    estimMSFE_H = zeros(tNum, hNum);

    parfor t = 1:tNum
        tLength = (t - 1) * tStep + tStart;
        estimMSFE_H_h = zeros(hNum, 1);
        for h = 1:hNum

            MA_H = params_H(h).MA(:);

            if ~p || (p == 1 && params_H(h).AR == 0)
                AR_H = [];
            else
                AR_H = params_H(h).AR(:);
            end

            p_H = length(AR_H);
            q_H = length(MA_H);
            
            % P_max_h is defined as P in point (ii)  of statement of Theorem 2.2 of Paper 1)
            r_star = max([p_H, q_H]);
            P_max_h = floor(tLength/h) + r_star;
            
            dBeta_reshaped_H = reshape(dBetaTotal(h).dBeta, p_H + q_H, p + q);
            V_beta_H = dBeta_reshaped_H * V_beta * dBeta_reshaped_H';
    
            [psiJac_H, piJac_H] = getJacobianPsiPi(AR_H, MA_H, params_H(h).Psi, params_H(h).Pi, P_max_h);

            for i = length(psiJac_H) + 1:P_max_h + 1
                psiJac_H(i, :) = 0;
            end
            for i = length(piJac_H) + 1:P_max_h + 1
                piJac_H(i, :) = 0;
            end
            % covMatrixPsiPi_H is the covariance matrix (see SIGMA_{XI_P} in point (ii) of statement of Theorem 2.2 of Paper 1)
            covMatrixPsiPi_H = [psiJac_H; piJac_H] * V_beta_H * [psiJac_H; piJac_H]';
            
            % estimation error for the given forecasting horizon h
            estimMSFE_H_h(h) = getEstimationErrorDaTa(params_H(h).Psi, params_H(h).sigma2, covMatrixPsiPi_H, 1, P_max_h, tLength, P_max_h);

            disp(h);
        end
        % estimation errors for all the hNum forecasting horizons 
        % and a given estimation sample length tLength
        estimMSFE_H(t, :) = estimMSFE_H_h';
    end
    
    disp('estimation error computation for the hybrid scheme is completed');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%                    COMPUTATION OF CHARACTERISTIC AND ESTIMATION ERRORS
    %%%%                    ASSOCIATED TO THE OPTIMAL HYBRID FORECASTING SCHEME
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % SELECTION OF THE BEST AGGREGATION SCHEME IN TERMS OF THE STEP/AGGREGATION
    % PERIOD VALUE IF H STEP HAS A FEW MULTIPLIERS, CALCULATION OF THE MINIMAL
    % ERROR VALUE AND OPTIMAL AGGREGATION PERIOD LENGTH

    % set the divisors of interest (see introduction of Section 3 of Paper 2)
    multArr = [2, 3, 4, 5, 6, 7, 8, 9, 10];
    multArrLength = length(multArr);
    estimMSFE_OH = zeros(tNum, 1);
    charactMSFE_OH = zeros(hNum, 1);
    curParam_H = struct('AR', [], 'MA', [], 'p', [], 'q', [], 'Psi', [], 'Pi', [], 'sigma2', []);
    MSFE_OH = struct('charactMSFE_OH', charactMSFE_OH, 'estimMSFE_OH', estimMSFE_OH, 'estimMSFE_OH_h', []);
    tempAggrStruct = struct('stepValues', [], 'MSFE_OH', MSFE_OH, 'curParam_H', curParam_H);

    for h = 1:hNum
        charactMSFE_OH = zeros(hNum, 1);
        curParam_H = struct('AR', [], 'MA', [], 'p', [], 'q', [], 'Psi', [], 'Pi', [], 'sigma2', []);
        MSFE_OH = struct('charactMSFE_OH', charactMSFE_OH, 'estimMSFE_OH', estimMSFE_OH, 'estimMSFE_OH_h', []);
        tempAggrStruct(h) = struct('stepValues', [], 'MSFE_OH', MSFE_OH, 'curParam_H', curParam_H);
        
        hLengthCur = hStart + (h - 1) * hStep;
        tempAggrStruct(h).MSFE_OH = [];
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
        
        % possible variants of aggregation periods/forecasting horizons
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
            tempAggrStruct(h).MSFE_OH(curStep).charactMSFE_OH = getCharForecastingErrorDaTa(params_H(curStep).Psi, params_H(curStep).sigma2, w);
            tempAggrStruct(h).curParam_H(curStep) = params_H(curStep);
            
            if (length(tempAggrStruct(h).curParam_H(curStep).AR) == 1 && tempAggrStruct(h).curParam_H(curStep).AR == 0)
                tempAggrStruct(h).curParam_H(curStep).p = 0;
                tempAggrStruct(h).curParam_H(curStep).AR = [];

            else
                tempAggrStruct(h).curParam_H(curStep).p = length(tempAggrStruct(h).curParam_H(curStep).AR);
            end
            
            tempAggrStruct(h).curParam_H(curStep).q = length(tempAggrStruct(h).curParam_H(curStep).MA);
        end
        
    end
    %
    % TOTAL ERROR ASSOCIATED TO THE OPTIMAL HYBRID FORECASTING SCHEME
    % 
    
    startMinValue = 10e6;
    
    totalOptError = zeros(tNum, hNum);
    optSteps = zeros(hNum, 1);
    charactMSFE_OH = zeros(tNum, hNum);
    for t = 1:tNum
        % We determine which aggregation period among the available
        % divisors is optimal for the given forecasting horizon of length h
        
        tLength = (t - 1) * tStep + tStart;
        
        for i = 1:hNum
            hLength = hStart + (i - 1) * hStep;
            curTempAggrStruct = tempAggrStruct(i);
            
            totalOptError(t, i) = startMinValue;
            
            for j = 1:length(curTempAggrStruct.stepValues)
                curStep = curTempAggrStruct.stepValues(j);
                
                hLengthCur = int8(hLength/curStep);
                curParam_OH = curTempAggrStruct.curParam_H(curStep);
                
                % P_max_h is defined as P in point (ii)  of statement of Theorem 2.2 of Paper 1)
                p_H = length(curParam_OH.AR);
                q_H = length(curParam_OH.MA);
                
                r_star = max([p_H, q_H]);
                
                w = zeros(1, hLengthCur);
                if strcmp(aggrType, 'stock')
                    w(1, end) = 1;
                else
                    w = ones(1, hLengthCur);
                end
                curdBeta_reshaped_OH = reshape(dBetaTotal(curStep).dBeta, curParam_OH.p + curParam_OH.q, p + q);
                curV_beta_OH = curdBeta_reshaped_OH * V_beta * curdBeta_reshaped_OH';
                
                P_max = floor(tLength/hLength) + r_star + hLengthCur - 1;
                P_max_h = floor(tLength/hLength) + r_star;% + hLengthCur - 1;
                
                [curPsiJac, curPiJac] = getJacobianPsiPi(curParam_OH.AR, curParam_OH.MA, curParam_OH.Psi, curParam_OH.Pi, P_max);
                
                for m = length(curPsiJac) + 1:P_max + 1
                    curPsiJac(m, :) = 0;
                end
                for m = length(curPiJac) + 1:P_max + 1
                    curPiJac(m, :) = 0;
                end
                curCovMatrixPsiPi_OH = [curPsiJac; curPiJac] * curV_beta_OH * [curPsiJac; curPiJac]';
                % compute the estimation error associated to the optimal
                % hybrid scheme with the aggregation period hLengthCur and
                % sample length floor(tLength/hLengthCur)
                tempAggrStruct(i).MSFE_OH(curStep).estimMSFE_OH_h = getEstimationErrorDaTa(curParam_OH.Psi, curParam_OH.sigma2, curCovMatrixPsiPi_OH, w, P_max_h, tLength, P_max);
                % compute the total forecasting error and choose the
                % optimal aggregation period for a given forecasting
                % horizon i
                if (totalOptError(t, i) > (tempAggrStruct(i).MSFE_OH(curStep).estimMSFE_OH_h + tempAggrStruct(i).MSFE_OH(curStep).charactMSFE_OH))
                   totalOptError(t, i) = tempAggrStruct(i).MSFE_OH(curStep).estimMSFE_OH_h  + tempAggrStruct(i).MSFE_OH(curStep).charactMSFE_OH;
                   optSteps(i) = curStep;
                end
            end
            charactMSFE_OH(t, i) = tempAggrStruct(i).MSFE_OH(optSteps(i)).charactMSFE_OH;
            
        end
    end
        
    
    % saving the results to the output variables
    charactMSFE_H = charactMSFE_H';

    totalMSFE_H = repmat(charactMSFE_H, tNum, 1) + estimMSFE_H;
    totalMSFE_AD = repmat(charactMSFE_AD, tNum, 1) + estimMSFE_AD;
    totalMSFE_OH = [];%totalOptError;
    
end

    
