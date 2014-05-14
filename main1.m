% this file is used to produce the results analogous to the ones given in
% Figure 1 of Paper 1
% evaluation of the total and characteristic errors associated to the following forecasting schemes:
% -All-disaggregated (AD)
% -All-aggregated (AA)
% -Hybrid (H)


clear;
matlabpool close force local;
matlabpool open;
% set up the ARMA coefficients
% UNCOMMENT THE ONES THAT ARE REQUIRED

% see THETA1 in Figure 1 of Paper 1
% AR = [];
% MA = [0, -0.54, 0, -0.24];

% see THETA2 in Figure 1 of Paper 1
AR = [];
MA = [-0.1, -0.54, 0.1, -0.24];


% get the ARMA order
p = length(AR);
q = length(MA);

% unconditional variance of the driving noise
sigma2 = 1;

% set the sample length for which the associated estimation errors have to be computed
% tStart fixes the smallest sample length considered
% if there is need to set the list of sample lengths for which the formulas
% have to be evaluated, then tStart fixes the lower bound of the interval,
% tEnd - higher bound and tStep fixes the step, respectively.
% Number of the sample lengths considered is computed via: tNum = floor((tEnd - tStart + tStep - hEnd) / tStep),
% where hEnd is the maximum forecasting horizon of interest.
% Then the i-th sample length is computed using: tLength(i) = (i - 1) * tStep + tStart,
% where i takes values from 1 up to tNum
tStart = 30;
tEnd = 5000;
tStep = 10;


% this parameter sets up the length of the Psi/Pi cut-off to be taken into
% account in the error computation
P_cut = tStart + max(p,q);


% setting the set of forecasting horizons of interest
% hStart fixes the smallest forecasting horizon considered
% if there is need to set the list of forecasting horizons for which the formulas
% have to be evaluated, then hStart fixes the lower bound of the interval,
% hEnd - higher bound and hStep fixes the step, respectively.
% Number of the forecasting horizons considered is computed via: hNum = (hEnd - hStart + hStep) / hStep.
% Then the i-th forecasting horizon is computed using: hLen(i) = hStart + (i - 1) * hStep,
% where i takes values from 1 up to hNum
%
% The exercise in Figure 1 of Paper 1 demonstrates the case when the
% aggregation period is equal to 2 and aggregation type is 'stock'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, we compute the total, estimation errors associated to the
% all-disaggregated forecasting scheme, hence we fix h=2 in the next
% commands as the forecasting horizon of interest
hStart_AD = 2;
hEnd_AD = 2;
hStep_AD = 1;

% set aggregation type
aggrType = 'stock';

% compute the total, characteristic and estimation errors, respectively,
% associated to the all-disaggregated forecasting scheme (AD)
[totalMSFE_AD, charactMSFE_AD, estimMSFE_AD] = getTheoreticalTotalErrorForARMA(AR, MA, sigma2, tStart, tEnd, tStep, ...
    hStart_AD, hEnd_AD, hStep_AD, P_cut);

% temporally aggregate the model ARMA(p,q) into ARMA(p,q*) with temporal
% aggregation period equal to 2
period_TA = 2;
% construct the aggregation vector w
w = zeros(1, period_TA);
w(end) = 1;

params_H.AR = AR;
params_H.MA = MA;
params_H.K = sigma2;

[AR_H, MA_H, sigma2_H] = temporalAggregationForARMA(p, q, params_H, w);
if ~(AR_H)
    AR_H = [];
end
if ~(MA_H)
    MA_H = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Second, we compute the total, estimation errors associated to the
% all-aggregated forecasting scheme, hence we fix h=1 in the next
% commands as the forecasting horizon of interest (since we are working at the aggregated level)
hStart_AA = 1;
hEnd_AA = 1;
hStep_AA = 1;

% compute the total, characteristic and estimation errors, respectively,
% associated to the all-aggregated forecasting scheme (AA)
% we divide the corresponding sample lengths by period_TA
[totalMSFE_AA, charactMSFE_AA, estimMSFE_AA] = getTheoreticalTotalErrorForARMA(AR_H, MA_H, sigma2_H, ...
    floor(tStart/period_TA), floor(tEnd/period_TA), floor(tStep/period_TA), ...
    hStart_AA, hEnd_AA, hStep_AA, P_cut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Third, we compute the total, estimation errors associated to the
% hybrid forecasting scheme, hence we fix h=1 in the next
% commands as the forecasting horizon of interest (since we are working at the aggregated level)
hStart_H = 1;
hEnd_H = 1;
hStep_H = 1;

% compute the total, characteristic and estimation errors, respectively,
% associated to the hybrid forecasting scheme (H)
[totalMSFE_H, charactMSFE_H, estimMSFE_H] = getTheoreticalTotalErrorForARMA(AR_H, MA_H, sigma2_H, ...
    tStart, tEnd, tStep, ...
    hStart_H, hEnd_H, hStep_H, P_cut);


% save results into the file, set up the filename
savefile = sprintf('Figure1_Paper1_ARMA%d%d_%s_example', p, q, aggrType);
save(savefile);