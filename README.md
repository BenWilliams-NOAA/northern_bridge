# northern_bridge

Convert GOA northern rockfish assessment from `ADMB` to `RTMB`

There are two primary folders, ADMB has the code and results for m22.1b, that was presented at the 2024 September groundfish Plan Team meeting. 
The R folder has a suite of scripts.

 - ADMB
 - R
   - bridge_data - same data used in m22.1b
   - bridge_pars - parameter outputs from m22.1b
   - bridge_model - two models `f `= m22.1b, `f1` adds priors on selectivity parameters
   - bridge - code for running the RTMB models
   - utils - self explanatory...


## Current status

The RTMB model `f` can effectively replicate ADMB `m22.1b` when parameter values are fixed.
However, when optimized the *M* goes to 0 and *q* is >5.
To address this *M* was fixed, and optimization was successful, however MCMC results were poor. 
Priors were added to selectivity (logistic) parameters using the starting values and CV=1 for lack of a better alternative. 
Results are similar to the ADMB model, but parameters do differ.
