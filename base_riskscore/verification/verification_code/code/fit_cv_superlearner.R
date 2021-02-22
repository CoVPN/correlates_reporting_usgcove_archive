#! /usr/bin/env/Rscript

# This script executes code under the heading of: Input to CV-SuperLearner 
# run on placebo arm data (CV-SL fit)

#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here("..", "..", ".."))
source(here::here("..", "..", "..", "_common.R"))
#-----------------------------------------------

# assume we will read in the iteration from the command line
args <- commandArgs(trailingOnly = TRUE)
this_seed <- as.numeric(args[[1]])
# library can be 'demo' or 'prod'
this_library <- args[[2]]
print(paste0("Running CV super learner for seed ", this_seed, " out of 10."))

# load random seeds
seeds <- readRDS(here::here("output", "seeds.rds"))
print(paste0("Actual seed set before run is ", seeds[this_seed],"."))

# From verification documents: input to CV.SuperLearner
# Y	outcome as vector derived from phase 1 dataset (ph1_placebo)	
# X	risk variables that have been scaled as data frame	
# family	binomial	
# obsWeights	All observations have equal weighting of 1 as this is phase 1 data	
# SL.library	SL library created as below based off combination of learners and screens	
# method	method.CC_nloglik	
# cvControl	list of 2 where first is V = V_outer, and other is stratifyCV which is TRUE	V_outer is outer validation fold equal to 5.
# innerCvControl	list of list containing V = total folds for inner validation	if np is less or equal to 30, then perform Leave-One-Out inner validation. Otherwise, perform 5-fold inner validation
# vimp	FALSE	This is not used as parameter in the CV-SL call but could be used in post-processing if TRUE.
library(SuperLearner)
# save outcome data
y <- readRDS(here::here("data_clean", "y.rds"))
x <- readRDS(here::here("data_clean", "x.rds"))

inner_validation_folds <- ifelse(sum(y) <= 30, sum(y) - 1, 5)
print(paste0("Super learners built using ", inner_validation_folds," folds of CV."))

# super learner library
if(this_library == "demo"){
  sl_library <- list(
    c("screen_all", "SL.mean"),
    c("screen_all", "SL.glm"),
    c("screen_glmnet", "SL.glm"),
    c("screen_univariate_logistic_pval", "SL.glm"),
    c("screen_highcor_random", "SL.glm")
  )
}else if(this_library == "prod"){
  sl_library <- list(
    c("screen_all", "SL.mean"),
    c("screen_all", "SL.glm"),
    c("screen_all", "SL.bayesglm"),
    c("screen_all", "SL.glm.interaction"),
    c("screen_glmnet", "SL.glm"),
    c("screen_univariate_logistic_pval", "SL.glm"),
    c("screen_highcor_random", "SL.glm"),
    c("screen_glmnet", "SL.bayesglm"),
    c("screen_univariate_logistic_pval", "SL.bayesglm"),
    c("screen_highcor_random", "SL.bayesglm"),
    c("screen_glmnet", "SL.glm.interaction"),
    c("screen_univariate_logistic_pval", "SL.glm.interaction"),
    c("screen_highcor_random", "SL.glm.interaction"),
    c("screen_glmnet", "SL.glmnet"),
    c("screen_univariate_logistic_pval", "SL.glmnet"),
    c("screen_highcor_random", "SL.glmnet"),
    c("screen_glmnet", "SL.gam"),
    c("screen_univariate_logistic_pval", "SL.gam"),
    c("screen_highcor_random", "SL.gam"),
    c("screen_glmnet", "SL.xgboost"),
    c("screen_univariate_logistic_pval", "SL.xgboost"),
    c("screen_highcor_random", "SL.xgboost"),
    c("screen_glmnet", "SL.cforest"),
    c("screen_univariate_logistic_pval", "SL.cforest"),
    c("screen_highcor_random", "SL.cforest")
  )
}

set.seed(seeds[this_seed])
cv_sl_fit <- CV.SuperLearner(
  Y = y,
  X = x,
  family = binomial(),
  obsWeights = rep(1, length(y)), 
  method = "method.CC_nloglik",
  cvControl = list(V = 5, stratifyCV = TRUE),
  innerCvControl = list(V = inner_validation_folds)
)

saveRDS(cv_sl_fit, file = here::here("output", paste0("cv_sl_fit_", this_seed, ".rds")))