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
source(here::here("code", "superlearner_library.R"))

inner_validation_folds <- ifelse(sum(y) <= 30, sum(y), 5)
print(paste0("Super learner built using ", inner_validation_folds," folds of CV."))

# NEED TO UPDATE BASED ON BHAVESH'S ADVICE
set.seed(this_seed)
system.time({
sl_fit <- SuperLearner(
  Y = y,
  X = data.frame(x),
  SL.library = sl_library,
  family = binomial(),
  method = "method.CC_nloglik",
  cvControl = list(V = inner_validation_folds,
                   stratifyCV = TRUE)
)
})

saveRDS(sl_fit, file = here::here("output", "sl_fit.rds"))