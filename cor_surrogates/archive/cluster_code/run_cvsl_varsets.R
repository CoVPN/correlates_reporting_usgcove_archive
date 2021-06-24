## load required libraries and functions
library(here)
library(methods)
library(SuperLearner)
library(e1071)
library(glmnet)
library(xgboost)
library(earth)
library(dplyr)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
#remotes::install_github("CoVPN/correlates_mockdata", auth_token = "")
library(COVIDcorr)
library(kyotil)
library(argparse)
library(vimp)
library(nloptr)
library(RhpcBLASctl)

# SL optimal surrogate analysis: Superlearner code requires computing environment with more than 10 cores!
num_cores <- parallel::detectCores()
if(num_cores < 10) stop("Number of cores on this computing environment are less than 10! Superlearner code needs atleast 11 cores to run smoothly.")

# Define code version to run
run_demo <- TRUE  # the demo version is simpler and runs faster!
run_prod <- FALSE # the prod version runs SL with the complete set of diverse learners and takes more time (around 5-7 hrs on a cluster machine!)

# source(paste0(here(), "/correlates_report/SL_estimated_optimal_surrogate/cluster_code/sl_screens.R")) # set up the screen/algorithm combinations
# source(paste0(here(), "/correlates_report/SL_estimated_optimal_surrogate/cluster_code/utils.R")) # get CV-AUC for all algs
source("sl_screens.R") # set up the screen/algorithm combinations
source("utils.R") # get CV-AUC for all algs

############ SETUP INPUT ####################### 
# Identify the endpoint variable
endpoint <- "EventIndPrimaryD57"
################################################    

# Create combined new dataset which has imputed values of demographics (for phase 1 data) from dat.covar.imp AND 
# imputed values for markers (for phase 2 data) from dat.wide.v
dat.ph1 <- COVIDcorr::dat.mock %>%
  filter(Perprotocol == 1) %>%
  filter(Trt == 1) %>% # consider only vaccine group
  mutate(Delta57overBbindSpike_2fold = ifelse(Day57bindSpike > (BbindSpike + log10(2)), 1, 0),
         Delta57overBbindSpike_4fold = ifelse(Day57bindSpike > (BbindSpike + log10(4)), 1, 0),
         Delta57overBbindRBD_2fold = ifelse(Day57bindRBD > (BbindRBD  + log10(2)), 1, 0),
         Delta57overBbindRBD_4fold = ifelse(Day57bindRBD > (BbindRBD  + log10(4)), 1, 0),
         Delta57overBpseudoneutid50_2fold = ifelse(Day57pseudoneutid50 > (Bpseudoneutid50  + log10(2)), 1, 0), 
         Delta57overBpseudoneutid50_4fold = ifelse(Day57pseudoneutid50 > (Bpseudoneutid50  + log10(4)), 1, 0), 
         Delta57overBpseudoneutid80_2fold = ifelse(Day57pseudoneutid80 > (Bpseudoneutid80  + log10(2)), 1, 0), 
         Delta57overBpseudoneutid80_4fold = ifelse(Day57pseudoneutid80 > (Bpseudoneutid80  + log10(4)), 1, 0), 
         Delta57overBliveneutid50_2fold = ifelse(Day57liveneutid50 > (Bliveneutid50  + log10(2)), 1, 0), 
         Delta57overBliveneutid50_4fold = ifelse(Day57liveneutid50 > (Bliveneutid50  + log10(4)), 1, 0), 
         Delta57overBliveneutid80_2fold = ifelse(Day57liveneutid80 > (Bliveneutid80  + log10(2)), 1, 0), 
         Delta57overBliveneutid80_4fold = ifelse(Day57liveneutid80 > (Bliveneutid80  + log10(4)), 1, 0)) 

dat.ph2 = dat.ph1 %>%
  filter(TwophasesampInd==1) %>%
  # Baseline Risk Factor includes only BRiskScore and Age (UPDATE THIS as mentioned in SAP)
  select(Ptid, Trt, BRiskScore, Age, EventIndPrimaryD57, wt,
         Day57bindSpike, Delta57overBbindSpike, Delta57overBbindSpike_2fold, Delta57overBbindSpike_4fold,
         Day57bindRBD, Delta57overBbindRBD, Delta57overBbindRBD_2fold, Delta57overBbindRBD_4fold,
         Day57pseudoneutid50, Delta57overBpseudoneutid50, Delta57overBpseudoneutid50_2fold, Delta57overBpseudoneutid50_4fold,
         Day57pseudoneutid80, Delta57overBpseudoneutid80, Delta57overBpseudoneutid80_2fold, Delta57overBpseudoneutid80_4fold,
         Day57liveneutid50, Delta57overBliveneutid50, Delta57overBliveneutid50_2fold, Delta57overBliveneutid50_4fold,
         Day57liveneutid80, Delta57overBliveneutid80, Delta57overBliveneutid80_2fold, Delta57overBliveneutid80_4fold) %>%
  filter(!is.na(Day57bindSpike), !is.na(Day57bindRBD), !is.na(Day57pseudoneutid50), !is.na(Day57pseudoneutid80), !is.na(Day57liveneutid50), !is.na(Day57liveneutid80))

Z_plus_weights <- dat.ph1 %>% 
  select(Ptid, EventIndPrimaryD57, wt, Trt, BRiskScore, Age) %>%
  #mutate(ptid = as.numeric(ptid)) %>% 
  filter(!is.na(EventIndPrimaryD57), !is.na(wt), !is.na(Trt), !is.na(BRiskScore), !is.na(Age))

###########################################################################
# Create combination scores across the 6 markers
dat.ph2 <- dat.ph2 %>% 
  left_join(get.pca.scores(dat.ph2 %>%
                 select(Ptid, Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80, Day57liveneutid50, Day57liveneutid80)), 
            by = "Ptid") %>%
  ## Similarly create nonlinear PCA columns and max signal diversity score
  mutate(FSDAM1 = sample(100, size = nrow(dat.ph2), replace = TRUE)/100,
         FSDAM2 = sample(100, size = nrow(dat.ph2), replace = TRUE)/100,
         max.signal.div.score = sample(1000, size = nrow(dat.ph2), replace = TRUE)/1000)

###########################################################################

# # Define treatment: either placebo/vaccine group for objective 3
# # Define endpoint: either y1/y2/y3 for objective 3
# treatment = 0
# endpoint = "y1"
# dat.ph1 <- dat.ph1 %>%  filter(trt == treatment)     # 784 records for placebo
# dat.ph2 <- dat.ph2 %>%  filter(trt == treatment)     # 155 records for placebo
# endpoint.ph1 = (dat.ph1 %>% select(matches(endpoint)))[[1]]
# endpoint.ph2 = (dat.ph2 %>% select(matches(endpoint)))[[1]]
# if(treatment == 0){
#   trt_string = "placebo"
# } else {
#   trt_string = "vaccine"
# }

###########################################################################

briskfactors <- dat.ph2 %>%
  select(BRiskScore, Age) %>%
  colnames()

markers <- dat.ph2 %>%
  select(Day57bindSpike:max.signal.div.score) %>%
  colnames()

#####################################################################################################################
## Create variable sets and set up X, Y for super learning
# Maternal enrollment variables are default in all sets

# 1. None (No markers; only maternal enrollment variables), phase 1 data
varset_baselineRiskFactors <- rep(FALSE, length(markers))

# 2-14
varset_bAbSpike <- create_varsets(markers, c("Day57bindSpike", "Delta57overBbindSpike", "Delta57overBbindSpike_2fold", "Delta57overBbindSpike_4fold"))
varset_bAbRBD <- create_varsets(markers, c("Day57bindRBD", "Delta57overBbindRBD", "Delta57overBbindRBD_2fold", "Delta57overBbindRBD_4fold"))
varset_pnabID50 <- create_varsets(markers, c("Day57pseudoneutid50", "Delta57overBpseudoneutid50", "Delta57overBpseudoneutid50_2fold", "Delta57overBpseudoneutid50_4fold"))
varset_pnabID80 <- create_varsets(markers, c("Day57pseudoneutid80", "Delta57overBpseudoneutid80", "Delta57overBpseudoneutid80_2fold", "Delta57overBpseudoneutid80_4fold"))
varset_lnabID50 <- create_varsets(markers, c("Day57liveneutid50", "Delta57overBliveneutid50", "Delta57overBliveneutid50_2fold", "Delta57overBliveneutid50_4fold"))
varset_lnabID80 <- create_varsets(markers, c("Day57liveneutid80", "Delta57overBliveneutid80", "Delta57overBliveneutid80_2fold", "Delta57overBliveneutid80_4fold"))
varset_bAb_pnabID50 <- create_varsets(markers, c("Day57bindSpike", "Delta57overBbindSpike", "Delta57overBbindSpike_2fold", "Delta57overBbindSpike_4fold",
                                                 "Day57bindRBD", "Delta57overBbindRBD", "Delta57overBbindRBD_2fold", "Delta57overBbindRBD_4fold",
                                                 "Day57pseudoneutid50", "Delta57overBpseudoneutid50", "Delta57overBpseudoneutid50_2fold", "Delta57overBpseudoneutid50_4fold"))
varset_bAb_pnabID80 <- create_varsets(markers, c("Day57bindSpike", "Delta57overBbindSpike", "Delta57overBbindSpike_2fold", "Delta57overBbindSpike_4fold",
                                                 "Day57bindRBD", "Delta57overBbindRBD", "Delta57overBbindRBD_2fold", "Delta57overBbindRBD_4fold",
                                                 "Day57pseudoneutid80", "Delta57overBpseudoneutid80", "Delta57overBpseudoneutid80_2fold", "Delta57overBpseudoneutid80_4fold"))
varset_bAb_lnabID50 <- create_varsets(markers, c("Day57bindSpike", "Delta57overBbindSpike", "Delta57overBbindSpike_2fold", "Delta57overBbindSpike_4fold",
                                                 "Day57bindRBD", "Delta57overBbindRBD", "Delta57overBbindRBD_2fold", "Delta57overBbindRBD_4fold",
                                                 "Day57liveneutid50", "Delta57overBliveneutid50", "Delta57overBliveneutid50_2fold", "Delta57overBliveneutid50_4fold"))
varset_bAb_lnabID80 <- create_varsets(markers, c("Day57bindSpike", "Delta57overBbindSpike", "Delta57overBbindSpike_2fold", "Delta57overBbindSpike_4fold",
                                                 "Day57bindRBD", "Delta57overBbindRBD", "Delta57overBbindRBD_2fold", "Delta57overBbindRBD_4fold",
                                                 "Day57liveneutid80", "Delta57overBliveneutid80", "Delta57overBliveneutid80_2fold", "Delta57overBliveneutid80_4fold"))
varset_bAb_combScores <- create_varsets(markers, c("Day57bindSpike", "Delta57overBbindSpike", "Delta57overBbindSpike_2fold", "Delta57overBbindSpike_4fold",
                                                   "Day57bindRBD", "Delta57overBbindRBD", "Delta57overBbindRBD_2fold", "Delta57overBbindRBD_4fold",
                                                   "PC1", "PC2", "FSDAM1", "FSDAM2", "max.signal.div.score"))
varset_allMarkers <- create_varsets(markers, markers[1:24])
varset_allMarkers_combScores <- create_varsets(markers, markers[1:29])


varset_names <- c("1_baselineRiskFactors", 
                  "2_varset_bAbSpike", "3_varset_bAbRBD", "4_varset_pnabID50", "5_varset_pnabID80",
                  "6_varset_lnabID50", "7_varset_lnabID80", "8_varset_bAb_pnabID50", "9_varset_bAb_pnabID80", 
                  "10_varset_bAb_lnabID50", "11_varset_bAb_lnabID80", "12_varset_bAb_combScores", "13_varset_allMarkers",
                  "14_varset_allMarkers_combScores")

## set up a matrix of all 
varset_matrix <- rbind(varset_baselineRiskFactors, 
                       varset_bAbSpike, varset_bAbRBD, varset_pnabID50, varset_pnabID80,
                       varset_lnabID50, varset_lnabID80, varset_bAb_pnabID50, varset_bAb_pnabID80,
                       varset_bAb_lnabID50, varset_bAb_lnabID80, varset_bAb_combScores, varset_allMarkers, 
                       varset_allMarkers_combScores)

job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
this_var_set <- varset_matrix[job_id, ]
cat("\n Running ", varset_names[job_id], "\n")

#############################################################################################################
X_covars2adjust_ph2 <- dat.ph2 %>% select(all_of(c(briskfactors, markers)))   

## scale covars2adjust to have mean 0, sd 1 for all vars
for (a in colnames(X_covars2adjust_ph2)) {
  X_covars2adjust_ph2[[a]] <- scale(X_covars2adjust_ph2[[a]],
                                    center = mean(X_covars2adjust_ph2[[a]], na.rm = T), 
                                    scale = sd(X_covars2adjust_ph2[[a]], na.rm = T))    
}

# X_covars2adjust_ph2 <- X_covars2adjust_ph2 %>% 
#   select_if(function(x) any(!is.na(x))) # Drop column if it has 0 variance, and returns all NAN's from scale function. 

##############################################################################################################
# select data based on job_id
# if(job_id == 1){
#   X_markers_varset <- X_covars2adjust_ph2[1:9]
#   Y = endpoint.ph2
#   weights = dat.ph2$wt
#   # X_markers_varset <- X_covars2adjust_ph1
#   # Y = endpoint.ph1
#   # weights = rep(1, length(Y))
#   sl_lib <- c("SL.mean", "SL.glm","SL.glm.interaction","SL.bayesglm", "SL.step", "SL.glmnet","SL.gam","SL.cforest","SL.xgboost")
#   #sl_lib <- SL_library
# } else {
  X_markers_varset <- bind_cols(X_covars2adjust_ph2[1:2], 
                                X_covars2adjust_ph2[3:length(X_covars2adjust_ph2)] %>% select(names(X_covars2adjust_ph2[3:length(X_covars2adjust_ph2)])[this_var_set])) %>%
              select_if(function(x) any(!is.na(x))) # Drop column if it has 0 variance, and returned all NAN's from scale function. 
  Y = dat.ph2$EventIndPrimaryD57
  weights = dat.ph2$wt
  sl_lib <- SL_library
# }

# C <- rep(1, length(Y))
# Z <- data.frame(Y = Y, X_markers_varset %>% select(age.at.trt, age.at.trt.cat, bmi, mhsmopr, m.ast, child5, season, smoker, daycare))

treatmentDAT <- dat.ph2 %>% select(Ptid, Trt, wt, EventIndPrimaryD57, all_of(c(briskfactors, markers))) %>%
  filter(Trt == 1) %>%
  select(-Trt)

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
if (sum(dat.ph2$EventIndPrimaryD57) <= 30){
  V_inner <- length(Y) - 1
} else if(sum(dat.ph2$EventIndPrimaryD57) > 30){
  V_inner <- 5
  #V_inner <- length(Y) - 1
}
  

## ---------------------------------------------------------------------------------
## run super learner, with leave-one-out cross-validation and all screens
## do 10 random starts, average over these
## use assay groups as screens
## ---------------------------------------------------------------------------------
## ensure reproducibility
set.seed(20201202)
# seeds <- round(runif(1, 1000, 10000)) # do only one seed as trial
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
                           # z_lib = c("SL.glm", "SL.bayesglm", "SL.step", "SL.gam","SL.cforest"), # new arguments
                           z_lib = "SL.glm",
                           scale = "logit", # new argument
                           vimp = FALSE,
                           mc.cores = num_cores
)

# saveRDS(fits, paste0("~/COVIDcorrSAP/results/SL_estimated_optimal_surrogate/sl_optsurrogate_cvaucs.rds"))

cvaucs <- list()
cvfits <- list()

for(i in 1:length(seeds)) {
  cvaucs[[i]] = fits[[i]]$cvaucs
  cvfits[[i]] = fits[[i]]$cvfits
}

saveRDS(cvaucs, paste0("~/COVIDcorrSAP/results/SL_estimated_optimal_surrogate/cvsl_optsurrogate_cvaucs_", "EventIndPrimaryD57_", "vacc_", varset_names[job_id], ".rds"))
save(cvfits, file = paste0("~/COVIDcorrSAP/results/SL_estimated_optimal_surrogate/cvsl_optsurrogate_cvfits_", "EventIndPrimaryD57_", "vacc_", varset_names[job_id], ".rda"))




# X_mat = X_markers_varset
# family = "binomial"
# Z = Z_treatmentDAT
# z_lib = c("SL.glm", "SL.bayesglm", "SL.step", "SL.gam","SL.cforest")
# z_lib = "SL.glm"
# obsWeights = weights
# all_weights = all_ipw_weights_treatment
# scale = "logit"
# method = "method.CC_nloglik"
# cvControl = list(V = V_outer, stratifyCV = TRUE)
# innerCvControl = list(list(V = V_inner))
# vimp = FALSE
# mc.cores = num_cores
