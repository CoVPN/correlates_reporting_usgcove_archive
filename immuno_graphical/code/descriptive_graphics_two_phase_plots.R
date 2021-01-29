#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(COVIDcorr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(GGally)
library(SWIM)
library(scales)
library(dummies)
library(gridExtra)
library(PResiduals)

# produce geom_statistics w/ resampling-based covariate-adjusted Spearman
source(here("code", "ggally_cor_resample.R"))
source(here("code", "covid_corr_plot_functions.R"))
source(here("code", "params.R"))

# load cleaned data
dat.long.twophase.sample <- readRDS(here("data_clean",
                                         "long_twophase_data.rds"))
dat.twophase.sample <- readRDS(here("data_clean", "twophase_data.rds"))


#-----------------------------------------------
# PAIR PLOTS
#-----------------------------------------------
# - The correlation of each pair of Day 57 and delta antibody marker readouts
#   are compared  within each of the baseline strata subgroups, treatment arm,
#   and baseline serostratum.
# - Pairs plots/scatterplots and baseline strata-adjusted Spearman rank
#   correlations are used.
# - Again, we made the correlation plots for each serostratum and all
#   serostrata combined
#-----------------------------------------------

for (tt in 1:4) {
  for (bserostatus in 0:1) {
    for (trt in c(0, 1)) {
      subdat <- dat.twophase.sample %>%
        dplyr::filter(Bserostatus == bserostatus & Trt == trt)

      covid_corr_pairplots(
        plot_dat = subdat,
        time = times[tt + 1],
        assays = assays,
        strata = "Bstratum",
        weight = "wt",
        plot_title = paste0(c("D29", "D57", "D29 Fold-rise over D1",
                              "D57 Fold-rise over D1")[tt],
                            " Ab markers: baseline ",
                            ifelse(bserostatus, "positive", "negative"), ", ",
                            c("placebo", "vaccine")[trt + 1], " arm"),
        column_labels = labels.axis[tt + 1, 1:4] %>% unlist,
        filename = paste0(save.results.to, "/pairs_", times[tt + 1],
                          "_Markers_", bstatus.labels.2[bserostatus + 1],
                          c("_placebo_arm", "_vaccine_arm")[trt + 1], "_",
                          study.name, ".png")
      )
    }
  }
}

for (tt in 1:4) {
  for (bserostatus in 0:1) {
    subdat <- dat.twophase.sample %>%
      dplyr::filter(Bserostatus == bserostatus)

    covid_corr_pairplots(
      plot_dat = subdat,
      time = times[tt + 1],
      assays = assays,
      strata = "Bstratum",
      weight = "wt",
      plot_title = paste0(c("D29", "D57", "D29 Fold-rise over D1",
                            "D57 Fold-rise over D1")[tt],
                          " Ab markers: baseline ",
                          ifelse(bserostatus, "positive", "negative"),
                          ", placebo + vaccine arm"),
      column_labels = labels.axis[tt + 1, 1:4] %>% unlist,
      filename = paste0(save.results.to, "/pairs_",
                        c("Day29", "Day57", "Delta29overB",
                          "Delta57overB")[tt], "_Markers_",
                        bstatus.labels.2[bserostatus+1], "_", study.name,
                        ".png")
    )
  }
}

#-----------------------------------------------
# The correlation of each pair of Day 1 antibody marker readouts are compared
# within each of the baseline demographic subgroups and baseline serostratum,
# pooling over the two treatment arms. Pairs plots/scatterplots and Spearman
# rank correlations are used, where the plots and tests pooling over the
# baseline demographic subgroups use the case-deleted data set.
#-----------------------------------------------

for (bserostatus in 0:1) {
  for (trt in 0:1) {
    subdat <- dat.twophase.sample %>%
      dplyr::filter(Bserostatus == bserostatus & as.numeric(Trt) == trt)

    covid_corr_pairplots(
      plot_dat = subdat,
      time = "B",
      assays = assays,
      strata = "Bstratum",
      weight = "wt",
      plot_title = paste0("D1 Ab markers: baseline ",
                          ifelse(bserostatus, "positive", "negative"), ", ",
                          c("placebo", "vaccine")[trt + 1], " arm"),
      column_labels = labels.axis[tt + 1, 1:4] %>% unlist,
      filename = paste0(save.results.to, "/pairs_baselineMarkers_",
                        bstatus.labels.2[bserostatus + 1], "_",
                        c("placebo", "vaccine")[trt + 1], "_arm_", study.name,
                        ".png")
    )
  }
}

## pairplots by baseline serostatus
for (bserostatus in 0:1) {
  subdat <- dat.twophase.sample %>%
    dplyr::filter(Bserostatus == bserostatus)

  covid_corr_pairplots(
    plot_dat = subdat,
     time = "B",
     assays = assays,
     strata = "Bstratum",
     weight = "wt",
     plot_title = paste0("D1 Ab markers: baseline ",
                         ifelse(bserostatus, "positive", "negative"),
                         ", vaccine + placebo arm"),
     column_labels = labels.axis[tt + 1, 1:4] %>% unlist,
     filename = paste0(save.results.to, "/pairs_baselineMarkers_",
                       bstatus.labels.2[bserostatus + 1], "_",
                       study.name, ".png")
  )
}


#-----------------------------------------------
# RCDF PLOTS
#-----------------------------------------------
# - Reverse empirical cdf (rcdf) plots for the Baseline, Day 57, and
#   Baseline-subtracted Day 57 assay readouts, stratified by treatment group
#   and baseline serostatus
# - We made four ggplot objects, each for one assay, and combine them with
#   ggarrange()
#-----------------------------------------------

for (tt in seq_along(times)) {
  covid_corr_rcdf_facets(
    plot_dat = dat.long.twophase.sample,
    x = times[tt],
    facet_by = "assay",
    color = "trt_bstatus_label",
    weight = "wt",
    panel_titles = labels.title2[tt,] %>% unlist,
    axis_titles = labels.axis[tt,] %>% unlist,
    filename = paste0(save.results.to, "/Marker_RCDF_", times[tt],
                      "_trt_both_bstatus_both_", study.name, ".png")
  )
}

#-----------------------------------------------
# RCDF plot of four day 29 and day 57 assay readouts in one plot, with the
# line-types  distinguishing the baseline serostatus
#-----------------------------------------------

for (day in c("29", "57")) {
  covid_corr_rcdf(
    plot_dat = subset(dat.long.twophase.sample, Trt = "Vaccine"),
    x = paste0("Day", day),
    color = "assay",
    lty = "Bserostatus",
    weight = "wt",
    xlab = paste0("D", day,
                  " Binding Ab (IU/ml) / Pseudovirus nAb ID50 or ID80"),
    filename = paste0(save.results.to, "/Marker_Rcdf_Day", day,
                      "_trt_vaccine_bstatus_both_", study.name, ".png")
  )
}

#-----------------------------------------------
# RCDF plot of four day 29 / day 57 fold-rise over baseline assay readouts in
# one plot, with the line-types distinguishing the baseline serostatus
#-----------------------------------------------

for (day in c("29", "57")) {
  covid_corr_rcdf(
    plot_dat = subset(dat.long.twophase.sample, Trt = "Vaccine"),
    x = paste0("Delta", day, "overB"),
    color = "assay",
    lty = "Bserostatus",
    weight = "wt",
    xlab = paste0("D", day, " Fold-rise over D1 Binding Ab (IU/ml) / Pseudovirus nAb ID50 or ID80"),
    filename = paste0(save.results.to, "/Marker_Rcdf_Delta", day,
                      "overB_trt_vaccine_bstatus_both_", study.name, ".png")
  )
}


#-----------------------------------------------
# RCDF plot of four day 57 assay readouts in one plot for only the vaccine
# recipients with baseline negative sero viral status
#-----------------------------------------------

for (day in c("29", "57")) {
  covid_corr_rcdf(
    plot_dat = subset(dat.long.twophase.sample, Trt = "Vaccine"),
    x = paste0("Day", day),
    color = "assay",
    lty = NULL,
    weight = "wt",
    xlab = paste0("D", day,
                  " Binding Ab (IU/ml) / Pseudovirus nAb ID50 or ID80"),
    filename = paste0(save.results.to, "/Marker_Rcdf_Day", day,
                      "_trt_vaccine_bstatus_Neg_", study.name, ".png")
  )
}

for (day in c("29", "57")) {
  covid_corr_rcdf(
    plot_dat = subset(dat.long.twophase.sample, Trt = "Vaccine"),
    x = paste0("Delta", day, "overB"),
    color = "assay",
    lty = NULL,
    weight = "wt",
    xlab = paste0("D", day,
                  " Fold-rise over D1 Binding Ab (IU/ml) / Pseudovirus nAb ID50 or ID80"),
    filename = paste0(save.results.to, "/Marker_Rcdf_Delta", day,
                      "overB_trt_vaccine_bstatus_Neg_", study.name, ".png")
  )
}


#-----------------------------------------------
# SCATTER PLOTS
#-----------------------------------------------
# - Scatter plots for correlation between Day29/Day57/Delta29overB/Delta57overB
#   and baseline, stratified by treatment, assay type and baseline serostatus.
# - We created four ggplot objects, each for one assay.
# - Each ggplot object is the scatter plot of Day 57 assay readout or
#   increase/fold-raise from baseline versus the baseline assay readouts,
#   stratified by the treatment group.
#-----------------------------------------------

for (tt in seq_along(times)) {
  for (bserostatus in 1:2) {
    for (trt in 1:2) {
      covid_corr_scatter_facets(
        plot_dat = subset(dat.long.twophase.sample,
                          as.numeric(Bserostatus) == bserostatus &
                          as.numeric(Trt) == trt)[, c(times, "assay", "Trt",
                                                  "Bstratum", "wt")],
        x = "B",
        y = times[tt],
        facet_by = "assay",
        strata = "Bstratum",
        weight = "wt",
        panel_titles = labels.assays.short %>% unlist,
        x_axis_titles = labels.axis[1, ] %>% unlist,
        y_axis_titles = labels.axis[tt, ] %>% unlist,
        filename = paste0(save.results.to, "/scatterplots_", times[tt], "vB_",
                          bstatus.labels.2[bserostatus],
                          c("_placebo_arm_", "_vaccine_arm_")[trt], study.name,
                          ".png"))
    }
  }
}


#-----------------------------------------------
# BOX PLOTS
#-----------------------------------------------
# - Box plots across treatment groups.
# - For the box plots, we made a ggplot object for every assay and use the
#   ggarrange() function to combine the resulted plots.
#-----------------------------------------------

#-----------------------------------------------
# boxplots of assay readouts at D1, D29 and D57, versus treatment groups for
# baseline negative subjects (main group of interest)
#-----------------------------------------------

for (bstatus in 1:2) {
  for (tt in seq_along(times)) {
    covid_corr_boxplot_facets(
      plot_dat = subset(dat.long.twophase.sample, 
                        Bserostatus == bstatus.labels[bstatus]),
      x = "Trt", 
      y = times[tt], 
      color = "Trt",
      facet_by = "assay", 
      plot_LLOQ = (tt <= 3), 
      LLOQ = LLOQ,
      legend = c("Placebo", "Vaccine"),
      axis_titles_y = labels.axis[tt, ] %>% unlist,
      panel_titles = labels.title2[tt, ] %>% unlist,
      filename = paste0(save.results.to, "/boxplots_", times[tt], "_x_trt_", bstatus.labels.2[bstatus],
                        "_", study.name, ".png"))
  }
}

#-----------------------------------------------
# - Box plots of the assay readouts versus baseline sero-status, stratified by
#   treatment groups
# - Make seperate plots for Placebo and Vaccine arms
#-----------------------------------------------

for (trt in 1:2) {
  for (tt in seq_along(times)) {
    covid_corr_boxplot_facets(
      plot_dat = subset(dat.long.twophase.sample, as.numeric(Trt) == trt),
      x = "Bserostatus",
      y = times[tt],
      color = "Bserostatus",
      facet_by = "assay",
      plot_LLOQ = (tt <= 3),
      LLOQ = LLOQ,
      legend = c("Baseline Negative", "Baseline Positive"),
      axis_titles_y = labels.axis[tt, ] %>% unlist,
      panel_titles = labels.title2[tt, ] %>% unlist,
      filename = paste0(save.results.to, "/boxplots_", times[tt],
                        "_x_bstatus_", c("placebo_arm_", "vaccine_arm_")[trt],
                        study.name, ".png")

  }
}
