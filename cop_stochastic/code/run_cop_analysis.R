#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

# load packages
library(here)
library(tidyverse)
library(conflicted)
library(hal9001)
library(sl3)
library(SuperLearner)
library(txshift)

# set options
sl3_debug_mode()
options(sl3.pcontinuous = 0)
conflict_prefer("filter", "dplyr")
plan(multicore, workers = 4L)

# run helper scripts
source(here("code", "params.R"))
source(here("code", "sl_lrnr_libs.R"))

# use faster SL library for testing?
if (run_forrest_run) {
  sl_learner_regression <- Lrnr_sl$new(
    learners = list(mean_lrnr, glm_lrnr),
    metalearner = logistic_metalearner
  )
}

# run analysis for each marker at each time
#future_lapply(markers, function(marker) {
lapply(markers, function(marker) {
  # get timepoint for marker
  this_time <- marker_to_time[[marker]]

  # load estimation-ready data for this timepoint
  data_est <- readRDS(
    here("data_clean", paste0("data_est_", this_time, ".rds"))
  )

  # estimate stochastic VE
  mcop_msm <- msm_vimshift(
    Y = data_est$outcome,
    A = as.numeric(scale(data_est[[marker]])),
    W = data_est[, covariates],
    C_samp = data_est$samp_ind,
    V = c("W", "Y"),
    delta_grid = seq(-0.5, 0.5, 0.1),
    msm_form = list(type = "linear", knot = NA),
    estimator = "tmle",
    weighting = "identity",
    ci_level = 0.95,
    ci_type = "marginal",
    # arguments passed to txshift()
    samp_fit_args = list(
      fit_type = "glm"
    ),
    g_exp_fit_args = list(
      fit_type = "sl",
      sl_learners_density = hose_hal_lrnr
    ),
    Q_fit_args = list(
      fit_type = "sl",
      #sl_learners = Lrnr_glm$new()
      sl_learners = sl_learner_regression
    ),
    eif_reg_type = "hal"
  )

  # TODO: save results? create plots?
})
#}, future.seed = TRUE)
