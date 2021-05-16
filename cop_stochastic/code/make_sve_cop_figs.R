#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

# load packages
library(here)
library(tidyverse)
library(latex2exp)
library(txshift)

# run helper scripts
source(here("code", "params.R"))

# load analysis results and make sumary figures for each marker at each time
lapply(markers, function(marker) {
  # get timepoint for marker and marker name
  this_time <- marker_to_time[[marker]]
  marker_name_short <- str_remove(marker, this_time)
  marker_name_long <- marker_to_name[[marker_name_short]]

  # load output data with risk-based and SVE results
  risk_results_msm <- readRDS(
    here("output", paste0("mcop_risk_", marker, ".rds"))
  )
  sve_results_msm <- readRDS(
    here("output", paste0("mcop_sve_", marker, ".rds"))
  )

  # plot for counterfactual infection risk in vaccinees
  p_risk_msm <- plot(risk_results_msm) +
    geom_hline(yintercept = risk_results_msm$param_est[delta == 0, psi],
               linetype = "dotted", colour = "red") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "red") +
    geom_errorbar(aes(ymin = ci_lwr, ymax = ci_upr), width = 0.2, size = 1,
                  linetype = "dashed") +
    geom_point(size = 14, alpha = 0.75) +
    coord_cartesian(ylim = c(0, NA)) +
    labs(
      x = paste("Shift in standardized", marker_name_long),
      y = "Stoch. Interv. Risk in Vaccinees",
      title = paste("Stoch. Interv. Risk of Symptomatic", "COVID-19 at Day",
                    parse_number(this_time)),
      subtitle = TeX(paste("working MSM summary:",
                           paste0("($\\hat{\\beta}_{TMLE}$ = ",
                           round(risk_results_msm$msm_est$param_est[2], 4),
                           ", ", "p-value = ",
                           round(risk_results_msm$msm_est$p_value[2], 4),
                           ")"))),
      caption = paste("
        Mean counterfactual COVID-19 infection risk across standardized shifts
        in", marker_name_long, "levels, summarized by projection of causal
        dose-response curve onto a linear working model.")
    ) +
    theme(
      legend.position = "none",
      text = element_text(size = 44),
      axis.text.x = element_text(colour = "black", size = 44, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 44)
    )
  ggsave(
    filename = here("figs", paste0("mcop_risk_", marker, ".pdf")),
    plot = p_risk_msm, height = 25, width = 28
  )

  # plot for counterfactual stochastic interventional vaccine efficacy
  p_sve_msm <- plot(sve_results_msm) +
    geom_hline(yintercept = sve_results_msm$param_est[delta == 0, psi],
               linetype = "dotted", colour = "red") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "red") +
    geom_errorbar(aes(ymin = ci_lwr, ymax = ci_upr), width = 0.2, size = 1,
                  linetype = "dashed") +
    geom_point(size = 14, alpha = 0.75) +
    coord_cartesian(ylim = c(NA, 1)) +
    labs(
      x = paste("Shift in standardized", marker_name_long),
      y = paste("Stoch. Interv. VE"),
      title = paste("Stoch. Interv. VE v.", "Symptomatic COVID-19 at Day",
                    parse_number(this_time)),
      subtitle = TeX(paste("working MSM summary:",
                           paste0("($\\hat{\\beta}_{TMLE}$ = ",
                           round(sve_results_msm$msm_est$param_est[2], 4),
                           ", ", "p-value = ",
                           round(sve_results_msm$msm_est$p_value[2], 4),
                           ")"))),
      caption = paste("
        Stochastic interventional vaccine efficacy v. COVID-19 infection across
        standardized shifts in", marker_name_long, "levels, summarized by
        projection of causal dose-response curve on a linear working model.")
    ) +
    theme(
      legend.position = "none",
      text = element_text(size = 44),
      axis.text.x = element_text(colour = "black", size = 44, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 44)
    )
  ggsave(
    filename = here("figs", paste0("mcop_sve_", marker, ".pdf")),
    plot = p_sve_msm, height = 25, width = 28
  )
})
