#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
source(here::here("..", "_common.R"))
#-----------------------------------------------

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
library(conflicted)
library(gam)
library(xgboost)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
load(here("output", "objects_for_running_SL.rda"))
source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
source(here("code", "utils.R")) # get CV-AUC for all algs

## solve cores issue
library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1)
print(blas_get_num_procs())
stopifnot(blas_get_num_procs() == 1)

## construct superlearner on placebo arm-----------------------
set.seed(20210216)
sl_riskscore_slfits <- SuperLearner(
  Y = Y, X = X_riskVars, family = "binomial",
  SL.library = SL_library, method = "method.CC_nloglik",
  cvControl = list(V = V_outer, stratifyCV = TRUE)
)

save(sl_riskscore_slfits, file = here("output", "sl_riskscore_slfits.rda"))

# Predict on vaccine arm
dat.ph1.vacc <- inputFile %>%
  filter(Perprotocol == 1 & Trt == 1) %>%
  # Keep only variables to be included in risk score analyses
  select(Ptid, Trt, all_of(endpoint), all_of(risk_vars)) %>%
  # Drop any observation with NA values in Ptid, Trt, or endpoint!
  drop_na(Ptid, Trt, all_of(endpoint))

X_covars2adjust_vacc <- dat.ph1.vacc %>%
  select(all_of(risk_vars))

# Impute missing values in any variable included in risk_vars using the mice package!
print("Make sure data is clean before conducting imputations!")
X_covars2adjust_vacc <- impute_missing_values(X_covars2adjust_vacc, risk_vars)

# Scale X_covars2adjust_vacc to have mean 0, sd 1 for all vars
for (a in colnames(X_covars2adjust_vacc)) {
  X_covars2adjust_vacc[[a]] <- scale(X_covars2adjust_vacc[[a]],
    center = mean(X_covars2adjust_vacc[[a]], na.rm = T),
    scale = sd(X_covars2adjust_vacc[[a]], na.rm = T)
  )
}

X_riskVars_vacc <- X_covars2adjust_vacc

pred_on_vaccine <- predict(sl_riskscore_slfits, newdata = X_riskVars_vacc, onlySL = TRUE)$pred %>%
  as.data.frame()

vacc <- bind_cols(
  dat.ph1.vacc %>% select(Ptid, all_of(endpoint)),
  pred_on_vaccine
) %>%
  rename(pred = V1) %>%
  # add AUC
  mutate(
    AUCchar = format(round(fast.auc(pred, EventIndPrimaryD57), 3), nsmall = 3),
    risk_score = log(pred / (1 - pred))
  )

write.csv(vacc, here("output", "vaccine_ptids_with_riskscores.csv"), row.names = FALSE)
