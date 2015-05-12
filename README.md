This repository contains functions that carry out the temporal aggregation of ARMA models and compute the forecasting error formulas explained in detail in the papers (referred as Paper 1 and Paper 2 in the code documentation)

Paper 1: Grigoryeva, L. and Ortega, J.-P. [2014] Hybrid forecasting with estimated temporally aggregated linear processes. To appear in the Journal of Forecasting. http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2148895

Paper 2: Grigoryeva, L. and Ortega, J.-P. [2014] Asymptotic forecasting error evaluation for estimated temporally aggregated linear processes. To appear in the International Journal of Computational Economics and Econometrics.

- main1.m allows to reproduce Figure 1 in Paper 1
- main2.m allows to reproduce the examples in Paper 2

REMARK: Notice that the parameters for solvers in functions temporalAggregation and temporalAggregationForARMA have to be tuned if the solver is displaying the error message. 
