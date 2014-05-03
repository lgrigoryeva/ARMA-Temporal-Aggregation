clear;

% set up the ARMA coefficients
AR = [0.8];    
MA = [-0.5, -0.54, 0.54, -0.24];    

% get the ARMA order
p = length(AR);
q = length(MA);

% conditional variance constant
sigma2 = 5;

% setting the length of the sample of a set of the lengths if 
% the result is needed for each of the lengths from the set
% tStep defines the step in the set
tStart = 15;
tEnd = 25;
tStep = 10;

% setting the set of forecasting horizons under question
% hStep defines the step in the set
hStart  = 1;
hEnd    = 10;
hStep   = 1;

% utility parameters for choosing the initial point for 
% nonlinear aggregation equations
pathLength = 10e4;
numTrials = 3;

% aggregation type
aggrType = 'stock';

% this parameter sets up the length of the Psi/Pi cut-of to be taken into
% account in the error computation
P_cut = tStart + max(p,q);

% computation of the total errors and characteristic errors associated to 
% the all-disaggregated, hybrid and optimal hybrid schemes
[totalMSFE_da, totalMSFE_ta, totalMSFE_ta_opt, ...
    charactMSFE_da, charactMSFE_ta, charactMSFE_ta_opt, optSteps] = getTheoreticalTotalError(AR, MA, sigma2, ...
    aggrType, tStart, tEnd, tStep, hStart, hEnd, hStep, pathLength, numTrials);

% save results into the file, set up the filename
savefile = sprintf('taResARMA%d%d_%s_example', p, q, aggrType);
save(savefile);