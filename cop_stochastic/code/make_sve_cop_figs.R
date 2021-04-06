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
library(conflicted)
conflict_prefer("filter", "dplyr")

# run helper scripts
source(here("code", "params.R"))

# load analysis results and make sumary figures for each marker at each time
lapply(markers, function(marker) {
  # load output data with risk-based and SVE results
  risk_results_msm <- readRDS(
    here("output", paste0("mcop_risk_", marker, ".rds"))
  )
  sve_results_msm <- readRDS(
    here("output", paste0("mcop_sve_", marker, ".rds"))
  )

  # TODO: make plots
  p_risk_msm <- plot(risk_results_msm) +
    geom_point(size = 6) +
    geom_hline(yintercept = risk_results_msm$param_est[delta == 0, psi],
               linetype = "dotted", colour = "red") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "blue") +
    labs(
      x = paste("Posited change in standardized RBD binding antibody",
                "(ID50)"),
      y = "Risk of symptomatic COVID-19 infection in vaccinees",
      title = "Estimates of mean symptomatic COVID-19 infection risks",
      subtitle = TeX(paste("working MSM summary:",
                           paste0("($\\hat{\\beta}_{TMLE}$ = ",
                           round(risk_results_msm$msm_est$param_est[2], 4),
                           ", ", "p-value = ",
                           round(risk_results_msm$msm_est$p_value[2], 4),
                           ")"))),
      caption = "
        Counterfactual COVID-19 infection risks across standardized shifts in
        RBD binding antibody activity levels. Projection of dose-response curve
        onto a linear working model provides evidence of decreased expected
        infection risk with increases in RBD binding antibody activity."
    ) +
    scale_x_continuous(breaks = seq(-1, 1, 0.1))
  p_risk_msm
})
