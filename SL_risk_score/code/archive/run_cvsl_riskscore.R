## load required libraries and functions
library("here")
library("methods")
library("SuperLearner")
library("e1071")
library("glmnet")
library("tidyverse")
## only run this if something has changed
# devtools::install_local("COVIDcorr.tar.gz")
library("COVIDcorr")
library("kyotil")
library("argparse")
library(vimp)
library(nloptr)
library(RhpcBLASctl)

num_cores <- parallel::detectCores()
#source(paste0(here(), "/correlates_report/SL_risk_score/cluster_code/sl_screens.R")) # set up the screen/algorithm combinations
#source(paste0(here(), "/correlates_report/SL_risk_score/cluster_code/utils.R")) # get CV-AUC for all algs

source("sl_screens.R") # set up the screen/algorithm combinations
source("utils.R") # get CV-AUC for all algs

# Consider only placebo data for risk score analysis
plac <- COVIDcorr::dat.mock %>%
  filter(Perprotocol == 1) %>%
  filter(Trt == 0) 

# Since in Moderna trial, the placebo arm has >400 cases or so, risk score will be built
# using not more than np/20 input variables, where np is the number of cases in the placebo arm. 
# Also, SL using 5-fold CV instead of leave one-out loss is advised. 

np <- floor(sum(plac$EventIndPrimaryD57)/20)

# Keep only variables to be included in risk score analyses
dat.ph1 <- plac %>%
  mutate(wt = 1,
         WhiteNonHispanic = ifelse(Hispanic==0 & Black==0 & Asian==0 &
                                     NatAmer==0 & PacIsl==0 & Multiracial==0 &
                                     Notreported==0 & Other==0 & Unknown==0, 1, 0)) %>%
  select(Ptid, Trt, EventIndPrimaryD57, wt, 
         MinorityInd, Hispanic, Black, Asian, NatAmer, PacIsl, Multiracial, 
         Notreported, Other, Unknown, WhiteNonHispanic,
         HighRiskInd, Sex, Age, BMI, BRiskScore) %>%
# Remove PacIsl as this race has fewer than 10 cases!
  select(-PacIsl)

Z_plus_weights <- dat.ph1 %>% 
  select(Ptid, Trt, EventIndPrimaryD57, wt, 
         MinorityInd, Hispanic, Black, Asian, NatAmer, Multiracial, 
         Notreported, Other, Unknown, WhiteNonHispanic,
         HighRiskInd, Sex, Age, BMI, BRiskScore) %>%
  filter(!is.na(EventIndPrimaryD57), !is.na(MinorityInd), !is.na(Hispanic), !is.na(Black),
         !is.na(Asian), !is.na(NatAmer), !is.na(Multiracial), 
         !is.na(Notreported), !is.na(Other), !is.na(Unknown), !is.na(WhiteNonHispanic),
         !is.na(HighRiskInd), !is.na(Sex), !is.na(Age), !is.na(BMI), !is.na(BRiskScore))

X_covars2adjust <-  dat.ph1 %>%
  select(MinorityInd, Hispanic, Black, Asian, NatAmer, Multiracial, 
         Notreported, Other, Unknown, WhiteNonHispanic, 
         HighRiskInd, Sex, Age, BMI, BRiskScore)
  
###########################################################################
# No missing value to impute!
###########################################################################
## scale covars2adjust to have mean 0, sd 1 for all vars
for (a in colnames(X_covars2adjust)) {
  X_covars2adjust[[a]] <- scale(X_covars2adjust[[a]],
                                center = mean(X_covars2adjust[[a]], na.rm = T),  
                                scale = sd(X_covars2adjust[[a]], na.rm = T))    
}


X_markers_varset <- X_covars2adjust
Y = plac$EventIndPrimaryD57
weights = dat.ph1$wt
sl_lib <- SL_library

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

## set up outer folds for cv variable importance; do stratified sampling
V_outer <- 5
V_inner <- 5

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
                           all_weights = all_ipw_weights_treatment,
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

saveRDS(cvaucs, paste0("~/COVIDcorrSAP/results/SL_risk_score/cvsl_riskscore_cvaucs.rds"))
save(cvfits, file = paste0("~/COVIDcorrSAP/results/SL_risk_score/cvsl_riskscore_cvfits.rda"))


