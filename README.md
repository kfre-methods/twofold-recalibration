# Two-fold recalibration of prediction models under selective predictor measurement

Generic R code for the test-based recalibration procedure described in the paper below.
When one predictor of a published model is measured in some patients only and imputed
otherwise, the procedure tests whether predictions with an observed and with an imputed
predictor require different corrections (pooled or two-fold recalibration), and selects
the most parsimonious correction by closed testing. It is implemented for time-to-event
models with competing risks, using pseudo-observations of the cumulative incidence.

`twofold_recalibration.R` contains the functions and an example on simulated data
(no real patient data). Required packages: survival, prodlim, geepack.

## Licence

Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0):
https://creativecommons.org/licenses/by-nc/4.0/

## Citation

If you use this code, please cite: Delord M, et al. Two-fold recalibration of prediction
models under selective predictor measurement. International Journal of Epidemiology (submitted).
