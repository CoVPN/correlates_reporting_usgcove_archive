#! /app/software/R/4.0.2-foss-2019b/bin/Rscript

# This script executes code under the heading of: Run Superlearner on placebo arm data (SL_fit) 

#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here("..", "..", ".."))
source(here::here("..", "..", "..", "_common.R"))
#-----------------------------------------------

# assume we will read in the iteration from the command line
args <- commandArgs(trailingOnly = TRUE)
this_seed <- 20210216
# library can be 'demo' or 'prod'
this_library <- args[[1]]
print(paste0("Running super learner for seed ", this_seed))

library(SuperLearner)
y <- readRDS(here::here("data_clean", "y.rds"))
x <- readRDS(here::here("data_clean", "x.rds"))

# super learner library
source(here::here("code", "sl_screen_fn.R"))
if(this_library == "demo"){
  sl_library <- list(
    c("SL.mean", "screen_all"),
    c("SL.glm", "screen_all"),
    c("SL.glm", "screen_glmnet"),
    c("SL.glm", "screen_univariate_logistic_pval"),
    c("SL.glm", "screen_highcor_random")
  )
}else if(this_library == "prod"){
  sl_library <- list(
    c("SL.mean", "screen_all"),
    c("SL.glm", "screen_all"),
    c("SL.bayesglm", "screen_all"),
    c("SL.glm.interaction", "screen_all"),
    c("SL.glm", "screen_glmnet"),
    c("SL.glm", "screen_univariate_logistic_pval"),
    c("SL.glm", "screen_highcor_random"),
    c("SL.bayesglm", "screen_glmnet"),
    c("SL.bayesglm", "screen_univariate_logistic_pval"),
    c("SL.bayesglm", "screen_highcor_random"),
    c("SL.glm.interaction", "screen_glmnet"),
    c("SL.glm.interaction", "screen_univariate_logistic_pval"),
    c("SL.glm.interaction", "screen_highcor_random"),
    c("SL.glmnet", "screen_glmnet"),
    c("SL.glmnet", "screen_univariate_logistic_pval"),
    c("SL.glmnet", "screen_highcor_random"),
    c("SL.gam", "screen_glmnet"),
    c("SL.gam", "screen_univariate_logistic_pval"),
    c("SL.gam", "screen_highcor_random"),
    c("SL.xgboost", "screen_glmnet"),
    c("SL.xgboost", "screen_univariate_logistic_pval"),
    c("SL.xgboost", "screen_highcor_random"),
    c("SL.cforest", "screen_glmnet"),
    c("SL.cforest", "screen_univariate_logistic_pval"),
    c("SL.cforest", "screen_highcor_random")
  )
}

inner_validation_folds <- ifelse(sum(y) <= 30, sum(y), 5)
print(paste0("Super learner built using ", inner_validation_folds," folds of CV."))

# NEED TO UPDATE BASED ON BHAVESH'S ADVICE
set.seed(this_seed)
sl_fit <- SuperLearner(
  Y = y,
  X = data.frame(x),
  SL.library = sl_library,
  family = binomial(),
  method = "method.CC_nloglik",
  cvControl = list(V = inner_validation_folds,
                   stratifyCV = TRUE)
)

saveRDS(sl_fit, file = here::here("output", "sl_fit.rds"))