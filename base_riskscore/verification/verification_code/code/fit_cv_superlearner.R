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

inner_validation_folds <- ifelse(sum(y) <= 30, sum(y), 5)
print(paste0("Super learners built using ", inner_validation_folds," folds of CV."))

# super learner library
source(here::here("code", "sl_screen_fn.R"))
source(here::here("code", "superlearner_library.R"))


set.seed(seeds[this_seed])
system.time({
cv_sl_fit <- CV.SuperLearner(
  Y = y,
  X = data.frame(x),
  SL.library = sl_library,
  family = binomial(),
  obsWeights = rep(1, length(y)), 
  method = "method.CC_nloglik",
  cvControl = list(V = 5, stratifyCV = TRUE),
  innerCvControl = list(list(V = inner_validation_folds))
)
})

saveRDS(cv_sl_fit, file = here::here("output", paste0("cv_sl_fit_", this_seed, ".rds")))