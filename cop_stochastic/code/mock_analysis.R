# packages, functions, options
library(here)
library(tidyverse)
library(latex2exp)
library(future)
library(sl3)
library(SuperLearner)
library(txshift)
source(here("..", "_common.R"))
source(here("code", "sl_libs.R"))
sl3_debug_mode()
options(sl3.pcontinuous = 0)
plan(transparent)

# load practice data and subset for second-stage sample
covid_ve <- read_csv(here("..", "data_raw", "ows_practice",
                          "covid_vetrial_mock.csv"))
colnames(covid_ve)[1] <- "ID"
covid_ve_tx <- covid_ve %>%
  dplyr::filter(Trt == 1, Perprotocol == 1, Bserostatus == 0)
head(covid_ve_tx)

#Mock data set attached. Some key variables:
#Trt = vaccine
#MinorityInd, HighRiskInd, Age = characteristics used for random subsample
#BRiskScore = another adjustment variable
#Day57liveneut = immune response of interest
#EventInd = outcome of interest
#So we’d restrict to Trt = 1
#A = Day57..live virus neutralizing antibody anti-Spike ID50
#W = MinorityInd, HighRiskInd, Age, BRiskScore
#And Y = EventInd
#Delta = !is.na(A)

#- Model sampling probs using saturated model of all components of W and Y.
#- Model other ones however you like (glm’s fine).
#- Key output is a figure like below that looks pretty.

# setup variables
ve_outcome <- covid_ve_tx$EventInd
ve_immuno_day57 <- covid_ve_tx$Day57liveneut
ve_immuno_day57_std <- as.numeric(scale(ve_immuno_day57))
ve_sample_ind <- as.numeric(!is.na(ve_immuno_day57))
ve_adjust_set <- covid_ve_tx[, c("MinorityInd", "Age", "BRiskScore",
                                 "HighRiskInd")]

# construct two-stage sampling model
sampling_data <- cbind(ve_sample_ind, ve_adjust_set)
glm_sampling <- glm(ve_sample_ind ~ .^4, data = sampling_data,
                    family = "binomial")
ipcw_fit_pred <- unname(predict(glm_sampling, type = "response"))
#ipcw_fit_pred[ipcw_fit_pred < 0.01] <- 0.01
ipc_weights <- ((ve_sample_ind == 1) / ipcw_fit_pred)[ve_sample_ind == 1]
ipcw_external <- list(pi_mech = ipcw_fit_pred, ipc_weights = ipc_weights)

# setup SL libraries g_A and Q_Y
sl_learner_regression <- Lrnr_sl$new(
  learners = list(mean_lrnr,
                  ranger_lrnr_base,
                  ranger_lrnr_ntrees50,
                  ranger_lrnr_ntrees100,
                  xgb_lrnr_base,
                  xgb50_lrnr,
                  xgb100_lrnr,
                  xgb300_lrnr,
                  bayesglm_sl_lrnr,
                  hal_lrnr_base,
                  hal_lrnr_custom,
                  glm_lrnr),
  metalearner = logistic_metalearner
)
sl_learner_density <- Lrnr_sl$new(
  learners = list(hose_glm_lrnr, hose_hal_lrnr),
  metalearner = Lrnr_solnp_density$new()
)

# run MSM analysis
ve_msm <- msm_vimshift(
  Y = ve_outcome,
  A = ve_immuno_day57_std,
  W = ve_adjust_set,
  C = ve_sample_ind,
  V = c("W", "Y"),
  delta_grid = seq(-0.9, 0.9, 0.1),
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
    #sl_learners_density = hese_glm_lrnr
    sl_learners_density = sl_learner_density
  ),
  Q_fit_args = list(
    fit_type = "sl",
    #sl_learners = Lrnr_glm$new()
    sl_learners = sl_learner_regression
  ),
  eif_reg_type = "hal",
  ipcw_fit_ext = ipcw_external
)
saveRDS(ve_msm, file = here("output", "covid_vetrial_mock_results.rds"))

# pretty plot for OWS presentation
p_ve_msm <- plot(ve_msm) +
  geom_point(size = 6) +
  geom_hline(yintercept = ve_msm$param_est[ve_msm$param_est$delta == 0]$psi,
             linetype = "dotted", colour = "red") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "blue") +
  labs(
    x = paste("Posited change in standardized live virus-neutralizing",
              "antibody (ID50)"),
    y = "Risk of symptomatic COVID-19 infection in vaccine recipients",
    title = paste("Estimates of mean symptomatic COVID-19 infection risks",
                  "under shifted live virus-neutralizing antibody"),
    subtitle = TeX(paste("with pointwise confidence intervals;",
                         "working marginal structural", "model summary",
                         paste0("($\\hat{\\beta}_{TMLE}$ = ",
                         round(ve_msm$msm_est$param_est[2], 4), ", ",
                         "p-value = ", round(ve_msm$msm_est$p_value[2], 3),
                         ")")))
  ) +
  annotate(geom = "text", size = 10, x = -0.4, y = 0.02,
           label = "Less potent antibodies") +
  annotate(geom = "text", size = 10, x = 0.4, y = 0.02,
           label = "More potent antibodies") +
  scale_y_continuous(breaks = seq(0, 0.025, by = 0.005)) +
  scale_x_continuous(breaks = seq(-0.9, 0.9, 0.1)) +
  coord_cartesian(xlim= c(-1, 1), ylim = c(0, 0.025))
p_ve_msm
ggsave2(filename = here("figs", "covid_mock_analysis.pdf"),
        plot = p_ve_msm)
