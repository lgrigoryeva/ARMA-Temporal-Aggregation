% this main file is used to launch the 
% evaluation of the total and characteristic errors associated to the following forecasting schemes:
% -All-disaggregated (AD)
% -Hybrid (H)
% -Optimal Hybrid (OH)
% The body contains a simple example 
% of how to use the code.


clear;
matlabpool close force local;
matlabpool open;
% set up the ARMA coefficients
AR = 0.8;    
MA = [-0.5, -0.54, 0.54, -0.24];    

% get the ARMA order
p = length(AR);
q = length(MA);

% unconditional variance of the driving noise
sigma2 = 5;

% set the sample length for which the associated estimation errors have to be computed
% tStart fixes the smallest sample length considered
% if there is need to set the list of sample lengths for which the formulas
% have to be evaluated, then tStart fixes the lower bound of the interval,
% tEnd - higher bound and tStep fixes the step, respectively.
% Number of the sample lengths considered is computed via: tNum = floor((tEnd - tStart + tStep - hEnd) / tStep),
% where hEnd is the maximum forecasting horizon of interest.
% Then the i-th sample length is computed using: tLength(i) = (i - 1) * tStep + tStart,
% where i takes values from 1 up to tNum

% In this case tNum = 1; tLength = 50;
tStart = 50;
tEnd = 70;
tStep = 20;

% setting the set of forecasting horizons of interest
% hStart fixes the smallest forecasting horizon considered
% if there is need to set the list of forecasting horizons for which the formulas
% have to be evaluated, then hStart fixes the lower bound of the interval,
% hEnd - higher bound and hStep fixes the step, respectively.
% Number of the forecasting horizons considered is computed via: hNum = (hEnd - hStart + hStep) / hStep.
% Then the i-th forecasting horizon is computed using: hLen(i) = hStart + (i - 1) * hStep,
% where i takes values from 1 up to hNum
hStart = 1;
hEnd = 10;
hStep = 1;

% utility parameters for choosing the initial point in the solution of the 
% nonlinear aggregation equations
pathLength = 10e4;%size of simulated path whose estimation provides initial conditions for the aggregated parameters
numTrials = 3;

% aggregation type
aggrType = 'stock';

% this parameter sets up the length of the Psi/Pi cut-off to be taken into
% account in the error computation
P_cut = tStart + max(p,q);

% computation of the total and characteristic errors associated to 
% the all-disaggregated (AD), hybrid (H) and optimal hybrid (OH) schemes
[totalMSFE_AD, totalMSFE_H, totalMSFE_OH, ...
    charactMSFE_AD, charactMSFE_H, charactMSFE_OH, optSteps, params_H] = getTheoreticalTotalError(AR, MA, sigma2, ...
    aggrType, tStart, tEnd, tStep, hStart, hEnd, hStep, pathLength, numTrials, P_cut);
% optSteps contains the aggregation periods corresponding to the OH scheme
% for each forecasting horizon
% params_H contains the parameters of the aggregated models


% save results into the file, set up the filename
savefile = sprintf('taResARMA%d%d_%s_example', p, q, aggrType);
save(savefile);