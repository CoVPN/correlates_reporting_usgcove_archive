# This script computes the cvauc of the fitted super learners

#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here("..", "..", ".."))
source(here::here("..", "..", "..", "_common.R"))
#-----------------------------------------------

library(vimp)
library(dplyr)

# Using SL fit object and scale = "identity", 
# get the CV-AUC estimate, SE and CI’s for each algorithm 
# (each individual learner, SL and discrete SL) as follows:
# For each of the 5 outer folds, run vimp::measure_auc  
# (vimp package version 2.1.10) to obtain cv-auc and eif 
# (estimated influence function)
# The cv-auc estimate is the average of cv-auc estimates 
# across the 5 outer folds.
# Get estimate of standard error using vimp::vimp_se function
# Get CI’s using vimp::vimp_ci function.

cv_auc_from_sl <- function(sl){
  V <- length(sl$folds)
  pt_est <- rep(NA, V)
  se_est <- rep(NA, V)
  for(v in 1:V){
    auc_v <- vimp::measure_auc(fitted_values = sl$SL.predict[sl$folds[[v]]],
                               y = sl$Y[sl$folds[[v]]])
    pt_est[v] <- auc_v$point_est
    se_est[v] <- as.numeric(vimp::vimp_se(auc_v))
  }
  cvauc_pt_est <- mean(pt_est)
  cvauc_se <- sqrt( mean(se_est^2) )
  ci <- vimp::vimp_ci(cvauc_pt_est, cvauc_se)

    	out <- cvAUC::ci.cvAUC(predictions = sl$SL.predict,
  	                       labels = sl$Y,
  	                       folds = sl$folds)
  	# replace ci with CI on logit scale
  	grad <- 1/(out$cvAUC - out$cvAUC^2)
    ci.logit <- plogis(qlogis(out$cvAUC) + sqrt(out$se^2 * grad^2) %o% c(-1.96, 1.96))
    out$ci <- as.numeric(ci.logit)
    return(out)
}

num_repeated_cv_seeds <- 1
rslt <- matrix(NA, nrow = num_repeated_cv_seeds, ncol = 5)
for(i in 1:num_repeated_cv_seeds){
  sl <- readRDS(here::here("output", paste0("cv_sl_fit_", i, ".rds")))
  rslt[i,] <- sl %>% cv_auc_from_sl %>% unlist
}

cv_auc_tab <- data.frame(rslt)[,1:4]
colnames(cv_auc_tab) <- c("cvAUC", "se", "ci.l", "ci.u")

saveRDS(cv_auc_tab, file = here::here("output", "cv_auc_tab.rds"))