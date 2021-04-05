# Receptor binding domain (RBD) binding antibody

# load results data from MSM analysis and plot
ve_msm <- readRDS(here("cop_stochastic", "output",
                       "ve_results_Day57bindRBD.rds"))
p_ve_msm <- plot(ve_msm) +
  geom_point(size = 6) +
  geom_hline(yintercept = ve_msm$param_est[ve_msm$param_est$delta == 0]$psi,
             linetype = "dotted", colour = "red") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "blue") +
  labs(
    x = paste("Posited change in standardized RBD binding antibody",
              "(ID50)"),
    y = "Risk of symptomatic COVID-19 infection in vaccinees",
    title = "Estimates of mean symptomatic COVID-19 infection risks",
    subtitle = TeX(paste("working MSM summary:",
                         paste0("($\\hat{\\beta}_{TMLE}$ = ",
                         round(ve_msm$msm_est$param_est[2], 4), ", ",
                         "p-value = ", round(ve_msm$msm_est$p_value[2], 4),
                         ")"))),
    caption = "
      Counterfactual COVID-19 infection risks across standardized shifts in RBD
      binding antibody activity levels. Projection of the dose-response curve
      onto a linear working model provides evidence of decreased expected
      infection risk with increases in RBD binding antibody activity."
  ) +
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.1))
p_ve_msm

# Spike protein antibody

# load results data from MSM analysis and plot
ve_msm <- readRDS(here("cop_stochastic", "output",
                       "ve_results_Day57bindSpike.rds"))
p_ve_msm <- plot(ve_msm) +
  geom_point(size = 6) +
  geom_hline(yintercept = ve_msm$param_est[ve_msm$param_est$delta == 0]$psi,
             linetype = "dotted", colour = "red") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "blue") +
  labs(
    x = paste("Posited change in standardized spike protein binding antibody",
              "(ID50)"),
    y = "Risk of symptomatic COVID-19 infection in vaccinees",
    title = "Estimates of mean symptomatic COVID-19 infection risks",
    subtitle = TeX(paste("working MSM summary:",
                         paste0("($\\hat{\\beta}_{TMLE}$ = ",
                         round(ve_msm$msm_est$param_est[2], 4), ", ",
                         "p-value = ", round(ve_msm$msm_est$p_value[2], 4),
                         ")"))),
    caption = "
     Counterfactual COVID-19 infection risks across standardized shifts
     in spike protein antibody activity levels. Projection of the dose-response
     curve onto a linear working model provides limited evidence of decreased
     expected infection risk with increases in spike protein antibody activity,
     though the trend flattens sharply beyond negative shifts."
  ) +
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.1))
p_ve_msm

# Live virus-neutralizing antibody

# load results data from MSM analysis and plot
ve_msm <- readRDS(here("cop_stochastic", "output",
                       "ve_results_Day57liveneutid80.rds"))
p_ve_msm <- plot(ve_msm) +
  geom_point(size = 6) +
  geom_hline(yintercept = ve_msm$param_est[ve_msm$param_est$delta == 0]$psi,
             linetype = "dotted", colour = "red") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "blue") +
  labs(
    x = paste("Posited change in standardized live virus-neutralizing",
              "antibody (ID50)"),
    y = "Risk of symptomatic COVID-19 infection in vaccinees",
    title = "Estimates of mean symptomatic COVID-19 infection risks",
    subtitle = TeX(paste("working MSM summary:",
                         paste0("($\\hat{\\beta}_{TMLE}$ = ",
                         round(ve_msm$msm_est$param_est[2], 4), ", ",
                         "p-value = ", round(ve_msm$msm_est$p_value[2], 4),
                         ")"))),
    caption = "
      Counterfactual COVID-19 infection risks across standardized shifts
      in live virus-neutralizing antibody activity levels. Projection of the
      dose-response curve onto a linear working model provides no evidence of
      of decreased expected infection risk with increases in live
      virus-neutralizing antibody activity."
  ) +
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.1))
p_ve_msm

# Pseudo-neutralizing antibody

# load results data from MSM analysis and plot
ve_msm <- readRDS(here("cop_stochastic", "output",
                       "ve_results_Day57pseudoneutid80.rds"))
p_ve_msm <- plot(ve_msm) +
  geom_point(size = 6) +
  geom_hline(yintercept = ve_msm$param_est[ve_msm$param_est$delta == 0]$psi,
             linetype = "dotted", colour = "red") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "blue") +
  labs(
    x = paste("Posited change in standardized pseudo-neutralizing antibody",
              "(ID50)"),
    y = "Risk of symptomatic COVID-19 infection in vaccinees",
    title = "Estimates of mean symptomatic COVID-19 infection risks",
    subtitle = TeX(paste("working MSM summary:",
                         paste0("($\\hat{\\beta}_{TMLE}$ = ",
                         round(ve_msm$msm_est$param_est[2], 4), ", ",
                         "p-value = ", round(ve_msm$msm_est$p_value[2], 4),
                         ")"))),
    caption = "
     Counterfactual COVID-19 infection risks across standardized shifts
     in pseudo-neutralizing antibody activity levels. Projection of the
     dose-response curve onto a linear working model provides limited evidence
     of decreased expected infection risk with increases in pseudo-neutralizing
     antibody activity."
  ) +
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.1))
p_ve_msm
