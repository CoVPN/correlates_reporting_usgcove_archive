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
library(conflicted)
conflict_prefer("filter", "dplyr")
library(mice)
############ Make sure latest COVIDcorr data package is installed!
latest_packageVersion <- sapply(strsplit(
  str_remove(list.files(path = paste0(here(), "/../R_packages"), pattern = ".tar.gz"), ".tar.gz"),
  "_"), `[`, 2)
if(latest_packageVersion != packageVersion("COVIDcorr"))
  install.packages(paste0(here(), "/R_packages/COVIDcorr_", latest_packageVersion, ".tar.gz"), repos = NULL, type ="source")
############ 
library(COVIDcorr)

# Risk score analysis: Superlearner code requires computing environment with more than 10 cores!
num_cores <- parallel::detectCores()
if(num_cores < 10) stop("Number of cores on this computing environment are less than 10! Superlearner code needs atleast 11 cores to run smoothly.")

# Define code version to run
run_demo <- FALSE  # the demo version is simpler and runs faster!
run_prod <- TRUE # the prod version runs SL with the complete set of diverse learners and takes more time (around 5-7 hrs on a cluster machine!)

#source(paste0(here(), "/correlates_report/SL_risk_score/cluster_code/sl_screens.R")) # set up the screen/algorithm combinations
#source(paste0(here(), "/correlates_report/SL_risk_score/cluster_code/utils.R")) # get CV-AUC for all algs
source("sl_screens.R") # set up the screen/algorithm combinations
source("utils.R") # get CV-AUC for all algs

############ SETUP INPUT ####################### 
# Identify the risk demographic variable names that will be used to compute the risk score
risk_vars <- c("MinorityInd", "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown", 
               "Black", "Asian", "NatAmer", "PacIsl", "WhiteNonHispanic", "Multiracial", 
               "Other", "Notreported", "Unknown", "HighRiskInd", "Sex", "Age", "BMI")

# Identify the endpoint variable
endpoint <- "EventIndPrimaryD57"
################################################    

# Consider only placebo data for risk score analysis
dat.ph1 <- COVIDcorr::dat.mock %>%
  filter(Perprotocol == 1 & Trt == 0) %>%
  mutate(wt = 1) %>% 
  # Keep only variables to be included in risk score analyses
  select(Ptid, Trt, all_of(endpoint), wt, all_of(risk_vars)) %>%  
  # Drop any observation with NA values for any of the selected variables
  drop_na()

# Limit total variables that will be included in models 
np <- sum(dat.ph1 %>% select(matches(endpoint)))
maxVar <- max(20, floor(np/20))

# Remove any risk_vars. For example, ones that have fewer than 10 cases (eg. EthnicityNotreported), or 
# ones with more than 5% (or 10%) missing values (if deemed unfit for imputation)! 
dat.ph1 <- drop_riskVars_with_fewer_cases(dat.ph1, risk_vars, endpoint)

# Update risk_vars
risk_vars <- dat.ph1 %>% select(-Ptid, -Trt, -all_of(endpoint), -wt) %>% colnames()

Z_plus_weights <- dat.ph1 %>% 
  select(Ptid, Trt, all_of(endpoint), wt, all_of(risk_vars)) %>%
  drop_na()

X_covars2adjust <-  dat.ph1 %>%
  select(all_of(risk_vars))
 
# Save ptids to merge with predictions later
risk_placebo_ptids <- dat.ph1 %>% select(Ptid, all_of(endpoint)) 

# Impute missing values in any variable included in risk_vars using the mice package!
X_covars2adjust <- impute_missing_values(X_covars2adjust, risk_vars)

# Scale covars2adjust to have mean 0, sd 1 for all vars
for (a in colnames(X_covars2adjust)) {
  X_covars2adjust[[a]] <- scale(X_covars2adjust[[a]],
                                center = mean(X_covars2adjust[[a]], na.rm = T),  
                                scale = sd(X_covars2adjust[[a]], na.rm = T))    
}

X_markers_varset <- X_covars2adjust
Y = dat.ph1 %>% pull(endpoint)
weights = dat.ph1$wt
sl_lib <- SL_library

# Create Z_treatmentDAT and C as inputs for vimp::measure_auc function
treatmentDAT <- dat.ph1 %>% select(-Trt)
# match the rows in treatmentDAT to get Z, C
all_cc_treatment <- Z_plus_weights %>%
  filter(Ptid %in% treatmentDAT$Ptid)
# pull out the participants who are NOT in the cc cohort and received the vaccine
all_non_cc_treatment <- Z_plus_weights %>%
  filter(!(Ptid %in% treatmentDAT$Ptid))
# put them back together
phase_1_data_treatmentDAT <- dplyr::bind_rows(all_cc_treatment, all_non_cc_treatment) %>%
  select(-Trt)
Z_treatmentDAT <- phase_1_data_treatmentDAT %>%
  select(-Ptid, -wt)
all_ipw_weights_treatment <- phase_1_data_treatmentDAT %>%
  pull(wt)
C <- (phase_1_data_treatmentDAT$Ptid %in% treatmentDAT$Ptid)


# set up outer folds for cv variable importance; do stratified sampling
V_outer <- 5
# set up inner folds based off number of cases 
if(np <= 30){
  V_inner <- length(Y) - 1
} else {
  V_inner <- 5  
}

## ensure reproducibility
set.seed(4747)
seeds <- round(runif(10, 1000, 10000)) # average over 10 random starts

##solve cores issue
library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1)
print(blas_get_num_procs())
stopifnot(blas_get_num_procs()==1)


fits <- parallel::mclapply(seeds, FUN = run_cv_sl_once, 
                           Y = Y, 
                           X_mat = X_markers_varset, 
                           family = "binomial",
                           obsWeights = weights,
                           all_weights = weights,
                           sl_lib = sl_lib,
                           method = "method.CC_nloglik",
                           cvControl = list(V = V_outer, stratifyCV = TRUE),
                           innerCvControl = list(list(V = V_inner)),
                           Z = Z_treatmentDAT, 
                           C = C, 
                           z_lib = "SL.glm",
                           scale = "logit", # new argument
                           vimp = FALSE,
                           mc.cores = num_cores
)

cvaucs <- list()
cvfits <- list()

for(i in 1:length(seeds)) {
  cvaucs[[i]] = fits[[i]]$cvaucs
  cvfits[[i]] = fits[[i]]$cvfits
}

saveRDS(cvaucs, "results/cvsl_riskscore_cvaucs.rds")
save(cvfits, file = "results/cvsl_riskscore_cvfits.rda")
save(risk_placebo_ptids, file = "results/risk_placebo_ptids.rda")
save(run_demo, run_prod, Y, X_markers_varset, weights, risk_vars, endpoint, maxVar, file = "results/objects_for_running_SL.rda")
