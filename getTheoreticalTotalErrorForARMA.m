function [totalMSFE, charactMSFE, estimMSFE] = getTheoreticalTotalErrorForARMA(AR, MA, sigma2, tStart, tEnd, tStep, ...
    hStart, hEnd, hStep, P_cut)
% This function evaluats the total, characteristic and estimation errors associated to the given forecasting scheme

% Output values: 
% totalMSFE - total forecasting error
% charactMSFE - characteristic forecasting error
% estimMSFE - estimation forecasting error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%                    INITIALIZATION BLOCK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
    % determine the order of ARMA
    p = length(AR);
    q = length(MA);
    r = max(p, q);
    
    %number of forecasting horizons
    hNum = (hEnd - hStart) / hStep + 1;

    % number of sample lengths of interest
    tNum = floor((tEnd - tStart) / tStep) + 1;

    % precisions and tolerance values
    precConvToMA    = 1.0e-12;
    recursStartLevel= 100;
    
    
    P_max = tEnd + hNum - 1 + r;
    if nargin < 10
        P_cut = 25;
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % CHARACTERISTIC FORECASTING ERROR
    
    charactMSFE = zeros(1, hNum);
    
    for h = 1:hNum
        hLen = hStart + (h - 1) * hStep;
        w = zeros(1, hLen);
        w(1, end) = 1;
        
        charactMSFE(h) = getCharForecastingErrorDaTa(Psi, sigma2, w);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %
    % ESTIMATION ERROR ASSOCIATED TO THE FORECASTING SCHEME
    % 
    
    estimMSFE = zeros(tNum, hNum);
        
    parfor t = 1:tNum
        estimMSFE_h = zeros(1, hNum);
        tLength     = (t - 1) * tStep + tStart;
        for h = 1:hNum
            hLen = hStart + (h - 1) * hStep;
            w = zeros(1, hLen);
            w(1, end) = 1;
            
            estimMSFE_h(h) = getEstimationErrorDaTa(Psi, sigma2, covMatrixPsiPi, w, min([P_cut, tLength]), tLength, P_max);
        end
        disp(t);
        estimMSFE(t, :) = estimMSFE_h;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %
    % TOTAL ERROR ASSOCIATED TO THE FORECASTING SCHEME
    % 
    totalMSFE = repmat(charactMSFE, tNum, 1) + estimMSFE;
end

    
