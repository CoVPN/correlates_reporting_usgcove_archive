# This script computes the cvauc of the fitted super learners

#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here("..", "..", ".."))
source(here::here("..", "..", "..", "_common.R"))
#-----------------------------------------------

library(cvAUC)
library(dplyr)

cv_auc_from_sl <- function(sl){
  	out <- cvAUC::ci.cvAUC(predictions = sl$SL.predict,
  	                       labels = sl$Y,
  	                       folds = sl$folds)
  	# replace ci with CI on logit scale
  	grad <- 1/(out$cvAUC - out$cvAUC^2)
    ci.logit <- plogis(qlogis(out$cvAUC) + sqrt(out$se^2 * grad^2) %o% c(-1.96, 1.96))
    out$ci <- as.numeric(ci.logit)
    return(out)
}

rslt <- matrix(NA, nrow = 10, ncol = 5)
for(i in 1:10){
  sl <- readRDS(here::here("output", paste0("cv_sl_fit_", i, ".rds")))
  rslt[i,] <- sl %>% cv_auc_from_sl %>% unlist
}

cv_auc_tab <- data.frame(rslt)[,1:4]
colnames(cv_auc_tab) <- c("cvAUC", "se", "ci.l", "ci.u")

saveRDS(cv_auc_tab, file = here::here("output", "cv_auc_tab.rds"))