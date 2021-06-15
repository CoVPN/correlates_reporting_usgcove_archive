###############################################################################
# Probably analyzing each of the four Day 57 antibody biomarkers (binding Ab to
# Spike protein, binding Ab to RBD, pseudo neut, live virus neut, the
# non-baseline subtracted variables), with one plot on results for each one.
###############################################################################

# packages, functions, options
library(here)
library(tidyverse)
library(latex2exp)
library(ggsci)
library(future)
library(hal9001)
library(sl3)
library(SuperLearner)
library(txshift)
source(here("..", "_common.R"))
source(here("code", "sl_libs.R"))
plan(multiprocess, workers = 4)
sl3_debug_mode()
options(sl3.pcontinuous = 0)

# longer names of immune markers of interest
marker_day57_longer_names <- list(
  Day57bindSpike = "spike protein binding antibody",
  Day57bindRBD = "RBD binding antibody",
  Day57pseudoneut = "pseudo-neutralizing antibody",
  Day57liveneut = "live virus-neutralizing antibody"
)

# load data and subset
covid_ve <- dat.mock %>%
  dplyr::filter(Trt == 1, Perprotocol == 1, Bserostatus == 0)
head(covid_ve)

# setup variables
ve_outcome <- covid_ve$EventIndPrimaryD57
ve_immuno_day57 <- covid_ve %>%
  select(
    "Day57bindSpike", "Day57bindRBD",
    "Day57pseudoneutid80", "Day57liveneutid80"
  )
ve_immuno_day57_std <- apply(ve_immuno_day57, 2, scale)
ve_sample_ind <- covid_ve$TwophasesampInd
apply(ve_immuno_day57, 2, function(immuno_day57) {
  stopifnot(all.equal(!is.na(immuno_day57), ve_sample_ind))
})
ve_adjust_set <- covid_ve %>%
  select(
    "MinorityInd", "Age", "BRiskScore", "HighRiskInd"
  )

# construct two-stage sampling model
sampling_data <- cbind(ve_sample_ind, ve_adjust_set)
glm_sampling <- glm(ve_sample_ind ~ .^4, data = sampling_data,
                    family = "binomial")
ipcw_fit_pred <- unname(predict(glm_sampling, type = "response"))
ipc_weights <- ((ve_sample_ind == 1) / ipcw_fit_pred)[ve_sample_ind == 1]
ipcw_external <- list(pi_mech = ipcw_fit_pred, ipc_weights = ipc_weights)

# run MSM analysis
lapply(seq_len(ncol(ve_immuno_day57_std)), function(col_idx) {
  # set immune biomarker for analysis
  marker_day57_name <- colnames(ve_immuno_day57_std)[col_idx]
  marker_day57_name_long <- marker_day57_longer_names[[marker_day57_name]]
  marker_day57_std <- ve_immuno_day57_std[, col_idx]

  # construct grid of shift estimates + summary MSM
  ve_msm <- msm_vimshift(
    Y = ve_outcome,
    A = marker_day57_std,
    W = ve_adjust_set,
    C = ve_sample_ind,
    V = c("W", "Y"),
    delta_grid = seq(-0.5, 0.5, 0.1),
    msm_form = list(type = "linear", knot = NA),
    estimator = "tmle",
    weighting = "identity",
    ci_level = 0.95,
    ci_type = "marginal",
    # arguments passed to txshift()
    ipcw_fit_args = list(
      fit_type = "external"
    ),
    g_fit_args = list(
      fit_type = "sl",
      sl_learners_density = hose_hal_lrnr
    ),
    Q_fit_args = list(
      fit_type = "sl",
      #sl_learners = Lrnr_glm$new()
      sl_learners = sl_learner_regression
    ),
    eif_reg_type = "hal",
    ipcw_fit_ext = ipcw_external
  )
  saveRDS(ve_msm,
          file = here("output",
                      paste0("ve_results_", marker_day57_name, ".rds")))

  # pretty plot for OWS presentation
  p_ve_msm <- plot(ve_msm) +
    geom_point(size = 6) +
    geom_hline(yintercept = ve_msm$param_est[ve_msm$param_est$delta == 0]$psi,
               linetype = "dotted", colour = "red") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "blue") +
    labs(
      x = paste("Posited change in standardized", marker_day57_name_long,
                "(ID50)"),
      y = "Risk of symptomatic COVID-19 infection in vaccine recipients",
      title = paste("Estimates of mean symptomatic COVID-19 infection risks",
                    "under shifted", marker_day57_name_long),
      subtitle = TeX(paste("with pointwise confidence intervals;",
                           "working marginal structural", "model summary",
                           paste0("($\\hat{\\beta}_{TMLE}$ = ",
                           round(ve_msm$msm_est$param_est[2], 4), ", ",
                           "p-value = ", round(ve_msm$msm_est$p_value[2], 4),
                           ")")))
    ) +
    #annotate(geom = "text", size = 10, x = -0.4, y = 0.02,
             #label = "Less potent antibodies") +
    #annotate(geom = "text", size = 10, x = 0.4, y = 0.02,
             #label = "More potent antibodies") +
    #scale_y_continuous(breaks = seq(0, 0.025, by = 0.005)) +
    scale_x_continuous(breaks = seq(-0.5, 0.5, 0.1))
    #coord_cartesian(xlim = c(-1, 1), ylim = c(0, 0.025))
  p_ve_msm
  ggsave2(filename = here("figs",
                          paste0("ve_summary_", marker_day57_name, ".pdf")),
         plot = p_ve_msm)
})
