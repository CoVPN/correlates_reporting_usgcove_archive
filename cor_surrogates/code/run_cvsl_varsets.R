#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
#-----------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
job_id <- as.numeric(args[1])
DAY <- as.character(args[2])

## load required libraries and functions
library(tidyverse)
library(quadprog)
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

# Assign weight and twophasesampind vars based off marker timepoint analyses
if(DAY %in% c("Day57", "Both")){
  WEIGHT = "wt.D57"
  TWOPHASESAMPIND = "TwophasesampIndD57"
}
if(DAY == "Day29"){
  WEIGHT = "wt.D29"
  TWOPHASESAMPIND = "TwophasesampIndD29"
}
# print(paste0("JOBID is ", args[1]))
# print(paste0("Day is ", args[2]))


# Define code version to run
# the demo version is simpler and runs faster!
# the production version runs SL with a diverse set of learners
run_prod <- !grepl("Mock", study_name)

# get utility files
source(here("code", "day57or29analyses.R")) # set up analyses for markers
source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
source(here("code", "utils.R")) # get CV-AUC for all algs

# SL optimal surrogate analysis: Superlearner code requires computing environment with more than 10 cores!
num_cores <- parallel::detectCores()
# if(num_cores < 10) stop("Number of cores on this computing environment are less than 10! Superlearner code needs atleast 11 cores to run smoothly.")

############ SETUP INPUT #######################
# Read data_clean
data_name_updated <- sub(".csv", "_with_riskscore.csv", data_name)
if (file.exists(here::here("..", "data_clean", data_name_updated))) {
  dat.mock <- read.csv(here::here("..", "data_clean", data_name_updated))
  data_name = data_name_updated
} else {
  dat.mock <- read.csv(here::here("..", "data_clean", data_name))
}

briskfactors <- c("risk_score", "HighRiskInd", "MinorityInd")

if(DAY == "Day57"){
  markerVars <- c("Day57bindSpike", "Delta57overBbindSpike", "Delta57overBbindSpike_2fold", "Delta57overBbindSpike_4fold",
                  "Day57bindRBD", "Delta57overBbindRBD", "Delta57overBbindRBD_2fold", "Delta57overBbindRBD_4fold",
                  "Day57pseudoneutid50", "Delta57overBpseudoneutid50", "Delta57overBpseudoneutid50_2fold", "Delta57overBpseudoneutid50_4fold",
                  "Day57pseudoneutid80", "Delta57overBpseudoneutid80", "Delta57overBpseudoneutid80_2fold", "Delta57overBpseudoneutid80_4fold")
}else if(DAY == "Day29"){
  markerVars <- c("Day29bindSpike", "Delta29overBbindSpike", "Delta29overBbindSpike_2fold", "Delta29overBbindSpike_4fold",
                  "Day29bindRBD", "Delta29overBbindRBD", "Delta29overBbindRBD_2fold", "Delta29overBbindRBD_4fold",
                  "Day29pseudoneutid50", "Delta29overBpseudoneutid50", "Delta29overBpseudoneutid50_2fold", "Delta29overBpseudoneutid50_4fold",
                  "Day29pseudoneutid80", "Delta29overBpseudoneutid80", "Delta29overBpseudoneutid80_2fold", "Delta29overBpseudoneutid80_4fold")
}else if(DAY == "Both"){
  markerVars <- c("Day29bindSpike", "Delta29overBbindSpike", "Delta29overBbindSpike_2fold", "Delta29overBbindSpike_4fold",
                  "Day29bindRBD", "Delta29overBbindRBD", "Delta29overBbindRBD_2fold", "Delta29overBbindRBD_4fold",
                  "Day29pseudoneutid50", "Delta29overBpseudoneutid50", "Delta29overBpseudoneutid50_2fold", "Delta29overBpseudoneutid50_4fold",
                  "Day29pseudoneutid80", "Delta29overBpseudoneutid80", "Delta29overBpseudoneutid80_2fold", "Delta29overBpseudoneutid80_4fold",
                  
                  "Day57bindSpike", "Delta57overBbindSpike", "Delta57overBbindSpike_2fold", "Delta57overBbindSpike_4fold",
                  "Day57bindRBD", "Delta57overBbindRBD", "Delta57overBbindRBD_2fold", "Delta57overBbindRBD_4fold",
                  "Day57pseudoneutid50", "Delta57overBpseudoneutid50", "Delta57overBpseudoneutid50_2fold", "Delta57overBpseudoneutid50_4fold",
                  "Day57pseudoneutid80", "Delta57overBpseudoneutid80", "Delta57overBpseudoneutid80_2fold", "Delta57overBpseudoneutid80_4fold")
}

# Identify the endpoint variable
endpoint <- "EventIndPrimaryD57"
################################################    

# Create combined new dataset which has imputed values of demographics (for phase 1 data) from dat.covar.imp AND 
# imputed values for markers (for phase 2 data) from dat.wide.v
dat.ph1 <- dat.mock %>%
  filter(Perprotocol == 1) %>%
  filter(Trt == 1) %>% # consider only vaccine group
  createBinaryVars(DAY) %>%
  # Drop any observation with NA values in Ptid, Trt, briskfactors, endpoint and wt.D57
  drop_na(Ptid, Trt, all_of(briskfactors), all_of(endpoint), all_of(WEIGHT)) %>%
  #dropNAforweightDay(DAY) %>%
  arrange(desc(get(endpoint)))

dat.ph2 <- dat.ph1 %>%
  filter(get(TWOPHASESAMPIND) == TRUE) %>%
  select(Ptid, Trt, all_of(briskfactors), all_of(endpoint), all_of(WEIGHT), any_of(markerVars)) %>%
  dropNAforDayMarker(DAY) %>%
  #drop_na(Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80) %>%
  arrange(desc(get(endpoint)))

# Limit total variables that will be included in models 
nv <- sum(dat.ph2 %>% select(matches(endpoint)))

# Remove any predictor variables that are indicator variables and have fewer than 10  0's or 1's
dat.ph2 <- drop_riskVars_with_fewer_0s_or_1s(dat.ph2, c(briskfactors, markerVars))

# Update predictor variables
pred_vars <- dat.ph2 %>%
  select(-Ptid, -Trt, -all_of(endpoint), -all_of(WEIGHT)) %>% 
  colnames()

# Remove any baseline risk factors with more than 5% missing values. Impute the missing
# values for other risk variables using mice package!
dat.ph2 <- drop_riskVars_with_high_total_missing_values(dat.ph2, briskfactors)

# Update risk_vars
pred_vars <- dat.ph2 %>%
  select(-Ptid, -Trt, -all_of(endpoint), -all_of(WEIGHT)) %>%
  colnames()

# Save ptids to merge with predictions later
ph2_vacc_ptids <- dat.ph2 %>% select(Ptid, all_of(endpoint), all_of(WEIGHT))

Z_plus_weights <- dat.ph1 %>% 
  select(Ptid, all_of(endpoint), all_of(WEIGHT), Trt, all_of(briskfactors)) %>%
  # Drop any observation with NA values in Ptid, Trt, briskfactors, endpoint or wt.D57
  drop_na(Ptid, Trt, all_of(briskfactors), all_of(endpoint), all_of(WEIGHT)) 
  
###########################################################################
# Create combination scores across the 5 markers
dat.ph2 <- dat.ph2 %>% 
  left_join(get.combination.scores(dat.ph2, DAY), by = "Ptid")

markers <- dat.ph2 %>%
  select(-Ptid, -Trt, -risk_score, -HighRiskInd, -MinorityInd, -EventIndPrimaryD57, -all_of(WEIGHT)) %>%
  colnames()

###########################################################################
## Create variable sets and set up X, Y for super learning
# Baseline risk variables are default in all sets

# 1. None (No markers; only baseline risk variables), phase 2 data
varset_baselineRiskFactors <- rep(FALSE, length(markers))

# 2-12
varset_bAbSpike <- create_varsets(markers, grep('bindSpike', markers, value=TRUE))
varset_bAbRBD <- create_varsets(markers, grep('bindRBD', markers, value=TRUE))
varset_pnabID50 <- create_varsets(markers, grep('pseudoneutid50', markers, value=TRUE))
varset_pnabID80 <- create_varsets(markers, grep('pseudoneutid80', markers, value=TRUE))
#varset_lnabMN50 <- create_varsets(markers, grep('liveneutmn50', markers, value=TRUE)) 
varset_bAb_pnabID50 <- create_varsets(markers, grep(paste(c('bindSpike', 'bindRBD', 'pseudoneutid50'), 
                                                          collapse="|"), markers, value=TRUE))
varset_bAb_pnabID80 <- create_varsets(markers, grep(paste(c('bindSpike', 'bindRBD', 'pseudoneutid80'), 
                                                          collapse="|"), markers, value=TRUE))
# varset_bAb_lnabMN50 <- create_varsets(markers, grep(paste(c('bindSpike', 'bindRBD', 'liveneutmn50'), 
#                                                           collapse="|"), markers, value=TRUE))
varset_bAb_combScores <- create_varsets(markers, 
                                        grep(paste(c('bindSpike', 'bindRBD', 'PC', 'nlPCA', 'max.signal.div.score'), 
                                                   collapse="|"), markers, value=TRUE))
varset_allMarkers <- create_varsets(markers, 
                                    grep(paste(c('bindSpike', 'bindRBD', 'pseudoneutid50', 'pseudoneutid80'), 
                                               collapse="|"), markers, value=TRUE))
varset_allMarkers_combScores <- create_varsets(markers, 
                                               grep(paste(c('bindSpike', 'bindRBD', 'pseudoneutid50', 'pseudoneutid80', 'PC', 'nlPCA', 'max.signal.div.score'), 
                                                          collapse="|"), markers, value=TRUE))

varset_names <- c("1_baselineRiskFactors",
                  "2_varset_bAbSpike", "3_varset_bAbRBD", "4_varset_pnabID50", "5_varset_pnabID80", 
                  "6_varset_bAb_pnabID50", "7_varset_bAb_pnabID80", 
                  "8_varset_bAb_combScores", "9_varset_allMarkers", "10_varset_allMarkers_combScores")

## set up a matrix of all
varset_matrix <- rbind(varset_baselineRiskFactors,
                       varset_bAbSpike, varset_bAbRBD, varset_pnabID50, varset_pnabID80, 
                       varset_bAb_pnabID50, varset_bAb_pnabID80, 
                       varset_bAb_combScores, varset_allMarkers, varset_allMarkers_combScores)

this_var_set <- varset_matrix[job_id, ]
cat("\n Running", varset_names[job_id], "variable set for", DAY, "markers! \n")

if(varset_names[job_id] == "1_noisyVariables"){
  X_covars2adjust_ph2 <- dat.ph2 %>% select(all_of(briskfactors), noiseVar1, noiseVar2, noiseVar3)
}else{
  X_covars2adjust_ph2 <- dat.ph2 %>% select(all_of(c(briskfactors, markers)))
}

## scale covars2adjust to have mean 0, sd 1 for all vars
for (a in colnames(X_covars2adjust_ph2)) {
  X_covars2adjust_ph2[[a]] <- scale(X_covars2adjust_ph2[[a]],
                                    center = mean(X_covars2adjust_ph2[[a]], na.rm = T),
                                    scale = sd(X_covars2adjust_ph2[[a]], na.rm = T))
}

markers_start <- length(briskfactors) + 1

if(job_id == 1){
  X_markers_varset <- X_covars2adjust_ph2[1:3] %>%
    select_if(function(x) any(!is.na(x))) # Drop column if it has 0 variance, and returned all NAN's from scale function.
}else{
  X_markers_varset <- bind_cols(X_covars2adjust_ph2[1:length(briskfactors)],
                                X_covars2adjust_ph2[markers_start:length(X_covars2adjust_ph2)][this_var_set]) %>%
    select_if(function(x) any(!is.na(x))) # Drop column if it has 0 variance, and returned all NAN's from scale function.
}

Y = dat.ph2 %>% pull(endpoint)
if(DAY %in% c("Day57", "Both")){
  weights = dat.ph2$wt.D57
}else if(DAY == "Day29"){
  weights = dat.ph2$wt.D29
}
  
sl_lib <- SL_library

if(varset_names[job_id] == "1_noisyVariables"){
  treatmentDAT <- dat.ph2 %>% select(Ptid, Trt, wt.D57, EventIndPrimaryD57, all_of(c(noiseVars, briskfactors, markers))) %>%
    filter(Trt == 1) %>%
    select(-Trt)
}else{
  treatmentDAT <- dat.ph2 %>% select(Ptid, Trt, all_of(WEIGHT), EventIndPrimaryD57, all_of(c(briskfactors, markers))) %>%
    filter(Trt == 1) %>%
    select(-Trt)
}

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
  select(-Ptid, -all_of(WEIGHT))
all_ipw_weights_treatment <- phase_1_data_treatmentDAT %>%
  pull(all_of(WEIGHT))
C <- (phase_1_data_treatmentDAT$Ptid %in% treatmentDAT$Ptid)

## set up outer folds for cv variable importance; do stratified sampling
V_outer <- 5
if (sum(dat.ph2 %>% pull(endpoint)) <= 25){
  V_inner <- length(Y) - 1
  maxVar <- 5
} else if(sum(dat.ph2 %>% pull(endpoint)) > 25){
  V_inner <- 5
  maxVar <- floor(nv/6)
}

## ---------------------------------------------------------------------------------
## run super learner, with leave-one-out cross-validation and all screens
## do 10 random starts, average over these
## use assay groups as screens
## ---------------------------------------------------------------------------------
## ensure reproducibility
set.seed(20210216)
seeds <- round(runif(10, 1000, 10000)) # average over 10 random starts

##solve cores issue
library(RhpcBLASctl)
#blas_get_num_procs()
blas_set_num_threads(1)
#print(blas_get_num_procs())
stopifnot(blas_get_num_procs()==1)

fits <- parallel::mclapply(seeds, FUN = run_cv_sl_once,
                           Y = Y,
                           X_mat = X_markers_varset,
                           family = "binomial",
                           obsWeights = weights,
                           all_weights = all_ipw_weights_treatment,
                           ipc_est_type = "ipw",
                           sl_lib = sl_lib,
                           method = "method.CC_nloglik",
                           cvControl = list(V = V_outer, stratifyCV = TRUE),
                           innerCvControl = list(list(V = V_inner)),
                           Z = Z_treatmentDAT,
                           C = C,
                           z_lib = "SL.glm",
                           scale = "identity", # new argument
                           vimp = FALSE,
                           mc.cores = num_cores
)


cvaucs <- list()
cvfits <- list()

for(i in 1:length(seeds)) {
  cvaucs[[i]] = fits[[i]]$cvaucs$aucs
  cvfits[[i]] = fits[[i]]$cvfits
}

saveRDS(cvaucs, file = here("output", paste0("CVSLaucs_vacc_", endpoint, "_", varset_names[job_id], "_", DAY, ".rds")))
save(cvfits, file = here("output", paste0("CVSLfits_vacc_", endpoint, "_", varset_names[job_id], "_", DAY, ".rda")))
save(ph2_vacc_ptids, file = here("output", paste0("ph2_vacc_ptids_", DAY, ".rda")))
save(run_prod, DAY, Y, dat.ph1, dat.ph2, weights, dat.mock, briskfactors, endpoint, maxVar,
     V_outer, file = here("output", paste0("objects_for_running_SL_", DAY, ".rda")))

