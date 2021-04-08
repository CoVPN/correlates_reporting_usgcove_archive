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
               linetype = "dashed", colour = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "blue") +
    geom_point(size = 4, alpha = 0.5) +
    labs(
      x = paste("Posited change in standardized", marker_name_long),
      y = "Risk of symptomatic COVID-19 infection in vaccinees",
      title = "Estimates of mean symptomatic COVID-19 infection risks",
      subtitle = TeX(paste("working MSM summary:",
                           paste0("($\\hat{\\beta}_{TMLE}$ = ",
                           round(risk_results_msm$msm_est$param_est[2], 4),
                           ", ", "p-value = ",
                           round(risk_results_msm$msm_est$p_value[2], 4),
                           ")"))),
      caption = paste("
        Counterfactual COVID-19 infection risks across standardized shifts
        in", marker_name_long, "levels, with projection of
        dose-response curve onto a linear working model for summarization.")
    ) +
    theme(
      legend.position = "bottom",
      legend.background =
        element_rect(fill = "gray90", size = 0.25, linetype = "dotted"),
      legend.title = element_blank(),
      text = element_text(size = 25),
      axis.text.x = element_text(colour = "black", size = 22, hjust = 1,
                                 angle = 30),
      axis.text.y = element_text(colour = "black", size = 22)
    )
  ggsave_custom(
    filename = here("figs", paste0("mcop_risk_", marker, ".pdf")),
    plot = p_risk_msm
  )

  # plot for counterfactual stochastic interventional vaccine efficacy
  p_sve_msm <- plot(sve_results_msm) +
    geom_hline(yintercept = sve_results_msm$param_est[delta == 0, psi],
               linetype = "dashed", colour = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "blue") +
    geom_point(size = 4, alpha = 0.5) +
    labs(
      x = paste("Posited change in standardized", marker_name_long),
      y = paste("Stochastic interventional vaccine efficacy against",
                "symptomatic COVID-19 infection"),
      title = paste("Estimates of mean vaccine efficacy against symptomatic",
                    "COVID-19 infection"),
      subtitle = TeX(paste("working MSM summary:",
                           paste0("($\\hat{\\beta}_{TMLE}$ = ",
                           round(sve_results_msm$msm_est$param_est[2], 4),
                           ", ", "p-value = ",
                           round(sve_results_msm$msm_est$p_value[2], 4),
                           ")"))),
      caption = paste("
        Stochastic counterfactual vaccine efficacy v. COVID-19 infection across
        standardized shifts in", marker_name_long, "levels, with projection
        of dose-response curve onto a linear working model for summarization.")
    ) +
    theme(
      legend.position = "bottom",
      legend.background =
        element_rect(fill = "gray90", size = 0.25, linetype = "dotted"),
      legend.title = element_blank(),
      text = element_text(size = 25),
      axis.text.x = element_text(colour = "black", size = 22, hjust = 1,
                                 angle = 30),
      axis.text.y = element_text(colour = "black", size = 22)
    )
  ggsave_custom(
    filename = here("figs", paste0("mcop_sve_", marker, ".pdf")),
    plot = p_sve_msm
  )
})
