This repository contains the R code used to develop and validate(internal and hold-out/on independent dataset) the RAISE (RA Infection Score Estimator) model for predicting the risk of serious infection (requiring hospitalisation) in rheumatoid arthritis patients initiating biologic or targeted synthetic DMARDs (b/tsDMARDs).
The script (prediction_model.R) includes:
Model derivation using a Cox proportional hazards model
Internal validation by bootstrapping (B = 500) and hold-out validation
Calibration assessment at 6, 12, 18, and 24 months
A prediction function (RAISE_predict()) to estimate infection risk for new patients
Individual-level SNDS data cannot be shared due to French data protection regulations.
