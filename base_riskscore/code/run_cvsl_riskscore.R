# Sys.setenv(TRIAL = "moderna_real")  
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
conflict_prefer("omp_set_num_threads", "RhpcBLASctl")
library(mice)
library(tidymodels)

# Define code version to run
# the demo version is simpler and runs faster!
# the production version runs SL with a diverse set of learners
run_prod <- !grepl("Mock", study_name)

# get utility files
source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
source(here("code", "utils.R")) # get CV-AUC for all algs

############ SETUP INPUT #######################
# Read in data file
inputFile <- read.csv(here::here("..", "data_clean", paste0(attr(config, "config"), "_data_processed.csv"))) 

# Identify the risk demographic variable names that will be used to compute the risk score
# Identify the endpoint variable
if(study_name_code == "COVE"){
  risk_vars <- c(
    "MinorityInd", "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown", 
    "Black", "Asian", "NatAmer", "PacIsl",  
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
    "URMforsubcohortsampling", "HighRiskInd", "HIVinfection", 
    "Sex", "Age", "BMI",
    "Country.X1", "Country.X2", "Country.X3", "Country.X4", "Country.X5", "Country.X6", "Country.X7", 
    "Region.X1", "Region.X2", 
    "CalDtEnrollIND.X1"
    )
  
  if(run_prod){
    risk_vars <- append(risk_vars, c("CalDtEnrollIND.X2", "CalDtEnrollIND.X3"))
  }
  
  endpoint <- "EventIndPrimaryIncludeNotMolecConfirmedD29"
  studyName_for_report <- "ENSEMBLE"

  # Create binary indicator variables for Country and Region
  inputFile <- inputFile %>%
    filter(!is.na(CalendarDateEnrollment)) %>%
    mutate(Sex.rand = sample(0:1, n(), replace = TRUE),
           Sex = ifelse(Sex %in% c(2, 3), Sex.rand, Sex),
           Country = as.factor(Country),
           Region = as.factor(Region),
           CalDtEnrollIND = case_when(CalendarDateEnrollment < 28 ~ 0,
                                      CalendarDateEnrollment >= 28 & CalendarDateEnrollment < 56 ~ 1,
                                      CalendarDateEnrollment >= 56 & CalendarDateEnrollment < 84 ~ 2,
                                      CalendarDateEnrollment >= 84 & CalendarDateEnrollment < 112 ~ 3,
                                      CalendarDateEnrollment >= 112 & CalendarDateEnrollment < 140 ~ 4,
                                      CalendarDateEnrollment >= 140 & CalendarDateEnrollment < 168 ~ 5),
           CalDtEnrollIND = as.factor(CalDtEnrollIND)) %>%
    select(-Sex.rand)
  
  rec <- recipe(~ Country + Region + CalDtEnrollIND, data = inputFile)
  dummies <- rec %>%
    step_dummy(Country, Region, CalDtEnrollIND) %>%
    prep(training = inputFile)
  inputFile <- inputFile %>% bind_cols(bake(dummies, new_data = NULL)) %>%
    select(-c(Country, Region, CalDtEnrollIND))
  names(inputFile)<-gsub("\\_",".",names(inputFile))
    
  # # Create interaction variables between Region and CalDtEnrollIND
  # rec <- recipe(EventIndPrimaryIncludeNotMolecConfirmedD29 ~., data = inputFile)
  # int_mod_1 <- rec %>%
  #   step_interact(terms = ~ starts_with("Region"):starts_with("CalDtEnrollIND"))
  # int_mod_1 <- prep(int_mod_1, training = inputFile)
  # inputFile <- bake(int_mod_1, inputFile)
  # names(inputFile)<-gsub("\\_",".",names(inputFile))
  # if(run_prod){
  #   risk_vars <- append(risk_vars, c("Region.X1.x.CalDtEnrollIND.X1", "Region.X1.x.CalDtEnrollIND.X2",
  #                                    "Region.X1.x.CalDtEnrollIND.X3",
  #                                    "Region.X2.x.CalDtEnrollIND.X1", "Region.X2.x.CalDtEnrollIND.X2",
  #                                    "Region.X2.x.CalDtEnrollIND.X3"))
  # }
}

################################################

# Consider only placebo data for risk score analysis
if("Riskscorecohortflag" %in% names(inputFile) & study_name_code != "COVE"){
  assertthat::assert_that(
    all(!is.na(inputFile$Riskscorecohortflag)), msg = "NA values present in Riskscorecohortflag in inputFile!"
  )
}else{
  if(study_name_code == "COVE"){
    inputFile <- inputFile %>%
      select(-Riskscorecohortflag) %>% # For Moderna, drop Riskscorecohortflag created in data_processing step and create a new simpler one!
      mutate(Riskscorecohortflag = ifelse(Bserostatus == 0 & Perprotocol == 1, 1, 0))
  }
  if(study_name_code == "ENSEMBLE"){
    inputFile <- inputFile %>%
      mutate(Riskscorecohortflag = ifelse(Bserostatus==0 & Perprotocol==1 & EarlyendpointD29==0 & EventTimePrimaryD29>=1, 1, 0)) 
  }
  
  assertthat::assert_that(
    all(!is.na(inputFile$Riskscorecohortflag)), msg = "NA values present in Riskscorecohortflag when created in inputFile!"
    )
}

dat.ph1 <- inputFile %>% filter(Riskscorecohortflag == 1 & Trt == 0)

dat.ph1 <- dat.ph1 %>%
  # Keep only variables to be included in risk score analyses
  select(Ptid, Trt, all_of(endpoint), all_of(risk_vars)) %>%
  # Drop any observation with NA values in Ptid, Trt, or endpoint!
  drop_na(Ptid, Trt, all_of(endpoint))

# Create table of cases in both arms (prior to Risk score analyses)
tab <- inputFile %>%
  filter(Riskscorecohortflag == 1) %>%
  drop_na(Ptid, Trt, all_of(endpoint)) %>%
  mutate(Trt = ifelse(Trt == 0, "Placebo", "Vaccine")) 

table(tab$Trt, tab %>% pull(endpoint)) %>%
  write.csv(file = here("output", "cases_prior_riskScoreAnalysis.csv"))
rm(tab)

# Derive maxVar: the maximum number of variables that will be allowed by SL screens in the models.
np <- sum(dat.ph1 %>% select(matches(endpoint)))
maxVar <- max(20, floor(np / 20))
all_risk_vars <- risk_vars

# Remove a variable if the number of cases in the variable = 1 subgroup is <= 3 or the number of cases in the variable = 0 subgroup is <= 3
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
save(run_prod, Y, X_riskVars, weights, inputFile, risk_vars, all_risk_vars, endpoint, maxVar,
     V_outer, studyName_for_report, file = here("output", "objects_for_running_SL.rda"))
