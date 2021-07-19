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
library(conflicted)
conflicted::conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
library(mice)

# Define code version to run
# the demo version is simpler and runs faster!
# the production version runs SL with a diverse set of learners
run_prod <- !grepl("Mock", study_name)

# get utility files
source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
source(here("code", "utils.R")) # get CV-AUC for all algs

############ SETUP INPUT #######################
# Read in data file
inputFile <- read.csv(here::here("..", "data_clean", data_name))

# Identify the risk demographic variable names that will be used to compute the risk score
# Identify the endpoint variable
if(study_name_code == "COVE"){
  risk_vars <- c(
    "MinorityInd", "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown", 
    "Black", "Asian", "NatAmer", "PacIsl", "WhiteNonHispanic", 
    "Multiracial", "Other", 
    "Notreported", "Unknown",
    "HighRiskInd", "Sex", "Age", "BMI"
  )
  
  endpoint <- "EventIndPrimaryD57"
  studyName_for_report <- "COVE"
}

if(study_name_code == "ENSEMBLE"){
  risk_vars <- c(
    "EthnicityHispanic","EthnicityNotreported", "EthnicityUnknown",
    "Black", "Asian", "NatAmer", "PacIsl", "Multiracial", "Notreported", "Unknown",
    "URMforsubcohortsampling", "HighRiskInd", "Sex", "Age", "BMI", "Country", 
    "HIVinfection", "CalendarDateEnrollment"
  )
  
  endpoint <- "EventIndPrimaryD29"
  studyName_for_report <- "ENSEMBLE"
}

################################################

# Consider only placebo data for risk score analysis
dat.ph1 <- inputFile %>%
  filter(Perprotocol == 1 & Trt == 0 & Bserostatus == 0) %>%
  # Keep only variables to be included in risk score analyses
  select(Ptid, Trt, all_of(endpoint), all_of(risk_vars)) %>%
  # Drop any observation with NA values in Ptid, Trt, or endpoint!
  drop_na(Ptid, Trt, all_of(endpoint))


# Derive maxVar: the maximum number of variables that will be allowed by SL screens in the models.
np <- sum(dat.ph1 %>% select(matches(endpoint)))
maxVar <- max(20, floor(np / 20))

# Remove any risk_vars that are indicator variables and have fewer than 10  0's or 1's
dat.ph1 <- drop_riskVars_with_fewer_0s_or_1s(dat = dat.ph1, 
                                             risk_vars = risk_vars,
                                             np = np)

# Update risk_vars
risk_vars <- dat.ph1 %>%
  select(-Ptid, -Trt, -all_of(endpoint)) %>%
  colnames()

# Remove any risk_vars with more than 5% missing values. Impute the missing
# values for other risk variables using mice package!
dat.ph1 <- drop_riskVars_with_high_total_missing_values(dat.ph1, risk_vars)

# Update risk_vars
risk_vars <- dat.ph1 %>%
  select(-Ptid, -Trt, -all_of(endpoint)) %>%
  colnames()

X_covars2adjust <- dat.ph1 %>%
  select(all_of(risk_vars))

# Save ptids to merge with predictions later
risk_placebo_ptids <- dat.ph1 %>% select(Ptid, all_of(endpoint))

# Impute missing values in any variable included in risk_vars using the mice package!
print("Make sure data is clean before conducting imputations!")
X_covars2adjust <- impute_missing_values(X_covars2adjust, risk_vars)

# # Check for missing values before and after imputation
# sapply(X_covars2adjust, function(x) sum(is.na(x)))

# Scale X_covars2adjust to have mean 0, sd 1 for all vars
for (a in colnames(X_covars2adjust)) {
  X_covars2adjust[[a]] <- scale(X_covars2adjust[[a]],
    center = mean(X_covars2adjust[[a]], na.rm = T),
    scale = sd(X_covars2adjust[[a]], na.rm = T)
  )
}

X_riskVars <- X_covars2adjust
Y <- dat.ph1 %>% pull(endpoint)

# set up outer folds for cv variable importance; do stratified sampling
V_outer <- 5

# set up inner folds based off number of cases
if (np <= 30) {
  V_inner <- length(Y) - 1
} else {
  V_inner <- 5
}

## solve cores issue
library(RhpcBLASctl)
#blas_get_num_procs()
blas_set_num_threads(1)
#print(blas_get_num_procs())
stopifnot(blas_get_num_procs() == 1)


# run super learner ensemble
fits <- run_cv_sl_once(
  seed = 20210216,
  Y = Y,
  X_mat = X_riskVars,
  family = "binomial",
  sl_lib = SL_library,
  method = "method.CC_nloglik",
  cvControl = list(V = V_outer, stratifyCV = TRUE),
  innerCvControl = list(list(V = V_inner)),
  scale = "identity",
  vimp = FALSE
)

cvaucs <- list()
cvaucs[[1]] <- fits$cvaucs
cvfits <- list()
cvfits[[1]] <- fits$cvfits

saveRDS(cvaucs, here("output", "cvsl_riskscore_cvaucs.rds"))
save(cvfits, file = here("output", "cvsl_riskscore_cvfits.rda"))
save(risk_placebo_ptids, file = here("output", "risk_placebo_ptids.rda"))
save(run_prod, Y, X_riskVars, weights, inputFile, risk_vars, endpoint, maxVar,
     V_outer, studyName_for_report, file = here("output", "objects_for_running_SL.rda"))
