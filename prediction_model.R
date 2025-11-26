###############################################################
# RAISE: Prediction model for serious infection risk in RA
# Author: Mounya Abboud
# Description: Development and internal validation of a Cox model 
# predicting serious infection (requiring hospitalisation) 
# in RA patients initiating b/tsDMARDs.
###############################################################

# Required packages
library(survival)
library(rms)
library(survAUC)
library(pec)
library(prodlim)
library(dplyr)

set.seed(123)

###############################################################
# 1. Model derivation (Cox proportional hazards)
###############################################################

# Define model formula
my_formula <- as.formula(
  Surv(time_to_event, infection) ~ 
    BEN_SEX_COD + age_cat + dose_category + molecule_of +
    infect_2ans_avant + pulmonary + diabetes + denutr +
    cancer + renal_failure + neuro + metho + 
    patho_card + MACE_history + htn + dyslipidemia +
    liver_pancreas + leflu + sulfa + hiv + obesity
)

# Fit the Cox model
cox_model <- coxph(my_formula, data = train_dataset)
summary(cox_model)

# Compute linear predictors (risk scores)
train_lp <- predict(cox_model, type = "lp")

# Discrimination (C-index)
train_concordance <- survConcordance(
  Surv(train_dataset$time_to_event, train_dataset$infection) ~ train_lp
)
train_concordance

###############################################################
# 2. Internal validation by bootstrapping (derivation cohort)
###############################################################

dd <- datadist(train_dataset)
options(datadist = "dd")

calib_model_train <- cph(
  my_formula,
  data = train_dataset,
  x = TRUE, y = TRUE, surv = TRUE
)

# Bootstrap validation (B = 500)
val <- validate(calib_model_train, method = "boot", B = 500, dxy = TRUE)
val  # optimism-corrected performance

###############################################################
# 3. Calibration at 6, 12, 18, and 24 months (derivation cohort)
###############################################################

cal_times <- c(6, 12, 18, 24)

calibrate_models_train <- lapply(cal_times, function(u) {
  calibrate(calib_model_train, method = "boot", B = 100, u = u)
})

# Example plot (24 months)
plot(calibrate_models_train[[4]],
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Predicted Probability",
     ylab = "Observed Probability",
     main = "Calibration - 24 months (Derivation Cohort)")

###############################################################
# 4. Validation on hold-out dataset (Validation Cohort)
###############################################################

test_lp <- predict(calib_model_train, newdata = test_dataset, type = "lp")

# C-index on validation cohort
test_concordance <- survConcordance(
  Surv(test_dataset$time_to_event, test_dataset$infection) ~ test_lp
)
test_concordance

# Calibration (same time points)
calibrate_models_test <- lapply(cal_times, function(u) {
  calibrate(fit = calib_model_train, method = "boot", u = u, B = 100, data = test_dataset)
})

# Example plot (24 months)
plot(calibrate_models_test[[4]],
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Predicted Probability",
     ylab = "Observed Probability",
     main = "Calibration - 24 months (Validation Cohort)")

###############################################################
# 5. Prediction function (example)
###############################################################
# Function to compute predicted risk of serious infection 
# for a new patient based on the final Cox model.
# Note: time_to_event is expressed in months in this dataset.
###############################################################
######################### my baseline survival hazard 
#Extract the baseline survival function
baseline_surv <- basehaz(calib_model_train, centered = TRUE)
baseline_surv

> Find the baseline hazard and calculate survival probability at 24 months
> baseline_hazard_24 <- baseline_surv[baseline_surv$time == 24, "hazard"]  
> baseline_survival_24 <- exp(-baseline_hazard_24)  
> baseline_survival_24
[1] 0.9400388

>Find the baseline hazard and calculate survival probability at 18 months
> baseline_hazard_18 <- baseline_surv[baseline_surv$time == 18, "hazard"]  
> baseline_survival_18 <- exp(-baseline_hazard_18)  
> baseline_survival_18
[1] 0.9527748

> Find the baseline hazard and calculate survival probability at 12 months
> baseline_hazard_12 <- baseline_surv[baseline_surv$time == 12, "hazard"]  # Hazard à 12 mois
> baseline_survival_12 <- exp(-baseline_hazard_12)  # Survie de base à 12 mois
> baseline_survival_12
[1] 0.9662098

> Find the baseline hazard and calculate survival probability at 6 months
> baseline_hazard_6 <- baseline_surv[baseline_surv$time == 6, "hazard"]  # Hazard à 6 mois
> baseline_survival_6 <- exp(-baseline_hazard_6)  # Survie de base à 6 mois
> baseline_survival_6
[1] 0.9811251
#######################


RAISE_predict <- function(newdata, time_months = 12) {
  # Compute the linear predictor for the new patient
  lp <- predict(cox_model, newdata = newdata, type = "lp")
  
  # Extract the baseline cumulative hazard from the fitted model
  base_surv <- basehaz(cox_model, centered = TRUE) #use the previous values for baseline hazard at 6,12,18,24 months
  
  # Identify the baseline cumulative hazard closest to the desired time (in months)
  surv_time <- base_surv$hazard[which.min(abs(base_surv$time - time_months))]
  
  # Compute predicted risk at the specified time horizon
  pred_risk <- 1 - exp(-exp(lp) * surv_time)
  
  # Return the predicted probability (numeric)
  return(pred_risk)
}


###############################################################
# Example: Predict 12-month infection risk for a new patient
###############################################################

new_patient <- data.frame(
  BEN_SEX_COD      = "F",            # Female
  age_cat          = "(58,67]",      # Age category
  dose_category    = ">=7.5 mg/day", # Corticosteroid dose category
  molecule_of      = "RITUXIMAB",    # Initiated b/tsDMARD
  infect_2ans_avant= 1,              # Prior serious infection (yes)
  pulmonary        = 1,              # Pulmonary disease (yes)
  diabetes         = 0,              # Diabetes (no)
  denutr           = 0,              # Severe malnutrition (no)
  cancer           = 0,              # Cancer (no)
  renal_failure    = 0,              # Renal failure (no)
  neuro            = 0,              # Neurologic disorder (no)
  metho            = 1,              # Methotrexate use (yes)
  patho_card       = 0,              # Non-atherothrombotic cardiac disease (no)
  MACE_history     = 0,              # Atherothrombotic disease (no)
  htn              = 1,              # Hypertension (yes)
  dyslipidemia     = 1,              # Dyslipidemia (yes)
  liver_pancreas   = 0,              # Liver/pancreas disease (no)
  leflu            = 0,              # Leflunomide (no)
  sulfa            = 0,              # Sulfasalazine (no)
  hiv              = 0,              # HIV infection (no)
  obesity          = 1               # Obesity (yes)
)

# Compute predicted risk at 12 months
RAISE_predict(new_patient, time_months = 12)

