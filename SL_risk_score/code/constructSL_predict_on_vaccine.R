# load required libraries and functions
library(tidyverse)
library(here)
library(methods)
library(SuperLearner)
library(e1071)
library(glmnet)
library(kyotil)
library(argparse)
library(vimp)
library(nloptr)
library(RhpcBLASctl)
library(aucm)
library(mice)
library(COVIDcorr)
library(conflicted)
conflict_prefer("filter", "dplyr")
load("results/objects_for_running_SL.rda")
source("sl_screens.R") # set up the screen/algorithm combinations
source("utils.R") 
sl_lib <- SL_library
## ensure reproducibility
set.seed(4747)

##solve cores issue
library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1)
print(blas_get_num_procs())
stopifnot(blas_get_num_procs()==1)

## construct superlearner on placebo arm-----------------------
sl_riskscore_slfits <- SuperLearner::SuperLearner(Y = Y, X = X_markers_varset, family = "binomial",
                           obsWeights = weights, SL.library = sl_lib,
                           method = "method.CC_nloglik")

save(sl_riskscore_slfits, file = "results/sl_riskscore_slfits.rda")

# Predict on vaccine arm
dat.ph1.vacc <- COVIDcorr::dat.mock %>%
  filter(Perprotocol == 1 & Trt == 1) %>%
  mutate(wt = 1) %>%
  # Keep only variables to be included in risk score analyses
  select(Ptid, Trt, all_of(endpoint), wt, all_of(risk_vars)) %>%
  # Drop any observation with NA values for any of the selected variables
  drop_na()

X_covars2adjust_vacc <-  dat.ph1.vacc %>%
  select(all_of(risk_vars))

# Impute missing values in any variable included in risk_vars using the mice package!
X_covars2adjust_vacc <- impute_missing_values(X_covars2adjust_vacc, risk_vars)

for (a in colnames(X_covars2adjust_vacc)) {
  X_covars2adjust_vacc[[a]] <- scale(X_covars2adjust_vacc[[a]],
                                center = mean(X_covars2adjust_vacc[[a]], na.rm = T),  
                                scale = sd(X_covars2adjust_vacc[[a]], na.rm = T))    
}

X_markers_varset_vacc <- X_covars2adjust_vacc

pred_on_vaccine = predict(sl_riskscore_slfits, newdata=X_markers_varset_vacc, onlySL = TRUE)$pred %>%
  as.data.frame()
vacc <- bind_cols(dat.ph1.vacc %>% select(Ptid, all_of(endpoint)),
                  pred_on_vaccine) %>%
  rename(pred = V1) %>%
  # add AUC
  mutate(AUCchar = format(round(fast.auc(pred, EventIndPrimaryD57),3), nsmall=3),
         risk_score = log(pred/(1-pred))) 

write.csv(vacc, "../output/vaccine_ptids_with_riskscores.csv", row.names = FALSE)                                  

