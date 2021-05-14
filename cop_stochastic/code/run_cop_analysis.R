#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

# load packages
library(here)
library(future)
library(future.apply)
library(progressr)
library(conflicted)
library(tidyverse)
library(sl3)
library(txshift)

# set options
sl3_debug_mode()
options(sl3.pcontinuous = 0)
conflict_prefer("filter", "dplyr")

# load data
data_name_amended <- paste0(str_remove(data_name, ".csv"),
                            "_with_riskscore.csv")
data_name_check <- file.exists(here("..", "data_clean", data_name_amended))

# run helper scripts
source(here("code", "params.R"))
source(here("code", "sl_lrnr_libs.R"))
source(here("code", "sve_utils.R"))

# parallelization and progress bar specification
plan(multicore, workers = availableCores())
handlers("txtprogressbar")

# add risk score to covariates list if defined
if (data_name_check) {
  covariates <- c(covariates, "risk_score")
}

# use faster SL library for testing
if (run_fast) {
  # simplify learner libraries for outcome and conditional density models
  sl_reg <- Lrnr_sl$new(
    learners = list(fglm, bayesglm, mean_y),
    metalearner = logistic_meta
  )
  sl_dens <- hose_glm
}

# run analysis for each marker at each time
with_progress({
  # instantiate progress bar
  p <- progressor(steps = length(markers))

  future_lapply(markers, function(marker) {
    # get timepoint for marker and marker name
    this_time <- marker_to_time[[marker]]
    marker_name <- str_remove(marker, this_time)

    # load estimation-ready data for this timepoint
    data_full <- readRDS(
      here("data_clean", paste0("data_est_", this_time, ".rds"))
    )
    data_est <- data_full %>%
      filter(Trt == 1) %>%
      select(-Trt)

    # estimate stochastic CoP based on risk
    mcop_risk_msm <- msm_vimshift(
      Y = data_est$outcome,
      A = as.numeric(scale(data_est[[marker]])),
      W = data_est[, covariates],
      C_samp = data_est$samp_ind,
      V = c("W", "Y"),
      delta_grid = delta_grid,
      msm_form = list(type = "linear", knot = NA),
      estimator = "tmle",
      weighting = "identity",
      ci_level = 0.95,
      ci_type = "marginal",
      # NOTE: arguments passed through to txshift()
      samp_fit_args = list(
        fit_type = "external"
      ),
      g_exp_fit_args = list(
        fit_type = "sl",
        sl_learners_density = sl_dens
      ),
      Q_fit_args = list(
        fit_type = "sl",
        sl_learners = sl_reg
      ),
      eif_reg_type = ifelse(run_fast, "glm", "hal"),
      # NOTE: need to pass in sampling probabilities, not weights
      samp_fit_ext = 1 / data_est$samp_wts
    )

    # transform risk estimates to the stochastic VE scale
    mcop_sve_msm <- sve_transform(
      mcop_risk_msm,
      data_full,
      weighting = "identity"
    )

    # increment the progress bar
    p(message = paste("Finished evaluating", marker_name, "at", this_time),
      class = "sticky")

    # save CoP risk results
    saveRDS(
      object = mcop_risk_msm,
      file = here("output", paste0("mcop_risk_", marker, ".rds"))
    )

    # save CoP SVE results
    saveRDS(
      object = mcop_sve_msm,
      file = here("output", paste0("mcop_sve_", marker, ".rds"))
    )
  }, future.seed = TRUE, future.conditions = "message")
})
