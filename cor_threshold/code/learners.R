#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------
source(here::here("code", "params.R"))
library(sl3)
library(SuperLearner)

bayesglm_sl_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.bayesglm", outcome_type = variable_type("binomial"))
# glm inluding two way interactions
lrnr_SL.inter <- make_learner(Lrnr_pkg_SuperLearner, "SL.glm.interaction", outcome_type = variable_type("binomial"))
# main term glm
lrnr_glm <- Lrnr_glm$new()
# main term lasso glm
lrnr_glmnet <- Lrnr_glmnet$new()
# Empirical mean estimator
lrnr_mean <- Lrnr_mean$new()
# baseline <- c("MinorityInd", "HighRiskInd", "Age", "BRiskScore", "age.geq.65")
interactions <- unlist(apply(combn(covariates, 2), 2, function(v) list(v)), recursive = F)
# Lasso with two way interactions
lrnr_glmnet_inter <- Pipeline$new(Lrnr_define_interactions$new(
  interactions
), lrnr_glmnet)
# bayes glm with interactions
lrnr_bayes_inter <- Pipeline$new(Lrnr_define_interactions$new(
  interactions
), bayesglm_sl_lrnr)


get_learner <- function(fast_analysis = T, include_interactions = F, covariate_adjusted = T) {
  if (!covariate_adjusted) {
    return(lrnr_mean)
  }
  if (fast_analysis) {
    if (include_interactions) {
      lrnr <- lrnr_glmnet_inter
    } else {
      lrnr <- lrnr_glmnet
    }
    return(lrnr)
  } else {
    if (include_interactions) {
      sl_list <- list(lrnr_bayes_inter, lrnr_glmnet, lrnr_mean, bayesglm_sl_lrnr, lrnr_glmnet_inter, lrnr_glm)
    } else {
      sl_list <- list(lrnr_glm, lrnr_glmnet, lrnr_mean, bayesglm_sl_lrnr)
    }
  }
  lrnr <- Lrnr_sl$new(sl_list, Lrnr_cv_selector$new(loss_loglik_binomial))
  return(lrnr)
}
