#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------
# install.packages(c("ggpubr", "GGally", "SWIM", "scales", "dummies",
# "gridExtra", "PResiduals"))
library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(GGally)
library(spatstat.geom)
library(scales)
library(dummies)
library(gridExtra)
library(PResiduals)

# produce geom_statistics w/ resampling-based covariate-adjusted Spearman
source(here("code", "ggally_cor_resample.R"))
source(here("code", "covid_corr_plot_functions.R"))
source(here("code", "params.R"))

set.seed(12345)
# load cleaned data
dat.long.twophase.sample <- readRDS(here(
  "data_clean",
  "long_twophase_data.rds"
))
dat.twophase.sample <- readRDS(here("data_clean", "twophase_data.rds"))

assay_lim <- readRDS(here("data_clean", "assay_lim.rds"))

tps <- c("Day29", "Day57", "Delta29overB", "Delta57overB")
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
print("Pair plots 1:")
for (tp in tps[tps %in% times]) {
  for (bserostatus in 0:1) {
    for (trt in c(0, 1)) {
      
      if(tp == "Day29"){
        tt <- 2
      }else if(tp == "Day57"){
        tt <- 3
      }else if(tp == "Delta29overB"){
        tt <- 4
      }else if(tp == "Delta57overB"){
        tt <- 5
      }
      
      subdat <- dat.twophase.sample %>%
        dplyr::filter(Bserostatus == bserostatus & Trt == trt)
      
      covid_corr_pairplots(
        plot_dat = subdat,
        time = tp,
        assays = assay_immuno,
        strata = "Bstratum",
        weight = "wt.subcohort",
        plot_title = paste0(
          c(
            "B", "D29", "D57", "D29 Fold-rise over D1",
            "D57 Fold-rise over D1"
          )[tt],
          " Ab markers: baseline ",
          ifelse(bserostatus, "positive", "negative"), ", ",
          c("placebo", "vaccine")[trt + 1], " arm"
        ),
        column_labels = labels.axis[tp, seq_along(assay_immuno)] %>% unlist(),
        height = max(1.3 * length(assay_immuno) + 0.1, 5.5),
        width = max(1.3 * length(assay_immuno), 5.5),
        filename = paste0(
          save.results.to, "/pairs_", tp,
          "_Markers_", bstatus.labels.2[bserostatus + 1],
          c("_placebo_arm", "_vaccine_arm")[trt + 1], "_",
          study_name, ".png"
        )
      )
    }
  }
}

#-----------------------------------------------
# The correlation of each pair of Day 1 antibody marker readouts are compared
# within each of the baseline demographic subgroups and baseline serostratum,
# pooling over the two treatment arms. Pairs plots/scatterplots and Spearman
# rank correlations are used, where the plots and tests pooling over the
# baseline demographic subgroups use the case-deleted data set.
#-----------------------------------------------
print("Pair plots 2:")
for (bserostatus in 0:1) {
  for (trt in 0:1) {
    subdat <- dat.twophase.sample %>%
      dplyr::filter(Bserostatus == bserostatus & as.numeric(Trt) == trt)
    
    covid_corr_pairplots(
      plot_dat = subdat,
      time = "B",
      assays = assay_immuno,
      strata = "Bstratum",
      weight = "wt.subcohort",
      plot_title = paste0(
        "D1 Ab markers: baseline ",
        ifelse(bserostatus, "positive", "negative"), ", ",
        c("placebo", "vaccine")[trt + 1], " arm"
      ),
      column_labels = labels.axis[1, seq_along(assay_immuno)] %>% unlist(),
      height = max(1.3 * length(assay_immuno) + 0.1, 5.5),
      width = max(1.3 * length(assay_immuno), 5.5),
      filename = paste0(
        save.results.to, "/pairs_baselineMarkers_",
        bstatus.labels.2[bserostatus + 1], "_",
        c("placebo", "vaccine")[trt + 1], "_arm_", study_name,
        ".png"
      )
    )
  }
}

print("Pair plots 3:")

## pairplots of assay readouts for multiple timepoints
## pairplots by baseline serostatus
for (bserostatus in 0:1) {
  for (trt in 0:1) {
    subdat <- dat.twophase.sample %>%
      dplyr::filter(Bserostatus == bserostatus & Trt == trt)
    
    for (aa in assay_immuno) {
      covid_corr_pairplots_by_time(
        plot_dat = subdat,
        times = intersect(c("B", "Day29", "Day57"), times),
        assay = aa,
        strata = "Bstratum",
        weight = "wt.subcohort",
        plot_title = paste0(
          labels.assays[aa], ": baseline ",
          ifelse(bserostatus, "positive ", "negative "),
          c("placebo", "vaccine")[trt + 1], " arm"
        ),
        column_labels = paste(c("D1", "D29", "D57")[c("B", "Day29", "Day57") %in% times], 
                              labels.axis[, aa][1]),
        filename = paste0(
          save.results.to, "/pairs_", aa, "_by_times_",
          bstatus.labels.2[bserostatus + 1], "_", c("placebo_", "vaccine_")[trt + 1],
          study_name, ".png"
        )
      )
    }
  }
}


#-----------------------------------------------
#-----------------------------------------------
# - Reverse empirical cdf (rcdf) plots for the Baseline, Day 57, and
#   Baseline-subtracted Day 57 assay readouts, stratified by treatment group
#   and baseline serostatus
# - We made four ggplot objects, each for one assay, and combine them with
#   ggarrange()
#-----------------------------------------------
print("RCDF 1:")
for (tp in tps[tps %in% times]) {
  covid_corr_rcdf_facets(
    plot_dat = dat.long.twophase.sample,
    x = tp,
    facet_by = "assay",
    color = "trt_bstatus_label",
    weight = "wt.subcohort",
    xlim = assay_lim[, tp, ],
    arrange_ncol = 3,
    arrange_nrow = ceiling(length(assay_immuno) / 3),
    panel_titles = labels.title2[tp, ] %>% unlist(),
    axis_titles = labels.axis[tp, ] %>% unlist(),
    filename = paste0(
      save.results.to, "/Marker_RCDF_", tp,
      "_trt_both_bstatus_both_", study_name, ".png"
    )
  )
}

#-----------------------------------------------
# RCDF plot of four day 29 and day 57 assay readouts in one plot, with the
# line-types  distinguishing the baseline serostatus
#-----------------------------------------------
dat.long.twophase.sample$assay_labels <-
  factor(dat.long.twophase.sample$assay,
         levels = assay_immuno,
         labels = labels.assays.short)
print("RCDF 2:")
# plot bAb and PsV assays separately
for (bAb in c(0, 1)) {
  
  if (bAb == 1) {
    rcdf_assays <- intersect(c("bindN", "bindSpike", "bindRBD"), assay_immuno) 
  } else {
    rcdf_assays <- assay_immuno[!assay_immuno %in% c("bindN", "bindSpike", "bindRBD")]
  }
  
  if (length(rcdf_assays) > 0) {
    
    for (tp in intersect(c("Day29", "Day57"), times)) {
      covid_corr_rcdf(
        plot_dat = subset(dat.long.twophase.sample, Trt == "Vaccine" & assay %in% rcdf_assays),
        x = tp,
        color = "assay_labels",
        lty = "Bserostatus",
        weight = "wt.subcohort",
        xlab = paste0(
          switch(tp, Day29 = "D29", Day57 = "D57"), " Ab Markers"
        ),
        xlim = c(min(assay_lim[assay_immuno %in% rcdf_assays, tp, 1]), 
                 max(assay_lim[assay_immuno %in% rcdf_assays, tp, 2])),
        xbreaks = seq(min(assay_lim[assay_immuno %in% rcdf_assays, tp, 1]), 
                      max(assay_lim[assay_immuno %in% rcdf_assays, tp, 2]), 
                      2),
        plot_title = paste0(switch(tp, Day29 = "Day 29", Day57 = "Day 57"), " Ab Markers"),
        filename = paste0(
          save.results.to, "/Marker_Rcdf_", c("nAb", "bAb")[bAb + 1], "_", tp,
          "_trt_vaccine_bstatus_both_", study_name, ".png"
        )
      )
    }
    
    #-----------------------------------------------
    # RCDF plot of four day 29 / day 57 fold-rise over baseline assay readouts in
    # one plot, with the line-types distinguishing the baseline serostatus
    #-----------------------------------------------
    print("RCDF 3:")
    for (tp in intersect(c("Delta29overB", "Delta57overB"), times)) {
      covid_corr_rcdf(
        plot_dat = subset(dat.long.twophase.sample, Trt == "Vaccine" & assay %in% rcdf_assays),
        x = tp,
        color = "assay_labels",
        lty = "Bserostatus",
        weight = "wt.subcohort",
        xlab = paste0(switch(tp, Delta29overB = "D29", Delta57overB = "D57"), 
                      " Fold-rise over D1 Ab Markers"),
        xlim = c(min(assay_lim[assay_immuno %in% rcdf_assays, tp, 1]), 
                 max(assay_lim[assay_immuno %in% rcdf_assays, tp, 2])),
        xbreaks = seq(min(assay_lim[assay_immuno %in% rcdf_assays, tp, 1]), 
                      max(assay_lim[assay_immuno %in% rcdf_assays, tp, 2]), 
                      2),
        plot_title = paste0(switch(tp, Delta29overB = "Day 29", Delta57overB = "Day 57"),
                            " over Baseline Ab Markers"),
        filename = paste0(
          save.results.to, "/Marker_Rcdf_", c("nAb", "bAb")[bAb + 1], "_", tp, 
          "_trt_vaccine_bstatus_both_", study_name, ".png"
        )
      )
    }
    
    
    #-----------------------------------------------
    # RCDF plot of four day 57 assay readouts in one plot for only the vaccine
    # recipients with baseline negative sero viral status
    #-----------------------------------------------
    print("RCDF 4:")
    for (bstatus in 1:2) {
      for (tp in intersect(c("Day29", "Day57"), times)) {
        covid_corr_rcdf(
          plot_dat = filter(dat.long.twophase.sample, Trt == "Vaccine", 
                            Bserostatus == bstatus.labels[bstatus],
                            assay %in% rcdf_assays),
          x = tp,
          color = "assay_labels",
          lty = NULL,
          weight = "wt.subcohort",
          xlab = paste0(
            switch(tp, Day29 = "D29", Day57 = "D57"), " Ab Markers"
          ),
          xlim = c(min(assay_lim[assay_immuno %in% rcdf_assays, tp, 1]), 
                   max(assay_lim[assay_immuno %in% rcdf_assays, tp, 2])),
          xbreaks = seq(min(assay_lim[assay_immuno %in% rcdf_assays, tp, 1]), 
                        max(assay_lim[assay_immuno %in% rcdf_assays, tp, 2]), 
                        2),
          plot_title = paste0(
            switch(tp, Day29 = "Day 29", Day57 = "Day 57"), " Ab Markers"
          ),
          filename = paste0(
            save.results.to, "/Marker_Rcdf_", c("nAb", "bAb")[bAb + 1], "_", tp,
            "_trt_vaccine_bstatus_", c("Neg", "Pos")[bstatus], "_", study_name, ".png"
          )
        )
      }
      
      for (tp in intersect(c("Delta29overB", "Delta57overB"), times)) {
        covid_corr_rcdf(
          plot_dat = filter(dat.long.twophase.sample, 
                            Trt == "Vaccine", 
                            Bserostatus == bstatus.labels[bstatus],
                            assay %in% rcdf_assays),
          x = tp,
          color = "assay_labels",
          lty = NULL,
          weight = "wt.subcohort",
          xlab = paste0(
            switch(tp, Delta29overB = "D29", Delta57overB = "D57"),
            " Fold-rise over D1 Ab Markers"
          ),
          xlim = c(min(assay_lim[assay_immuno %in% rcdf_assays, tp, 1]), 
                   max(assay_lim[assay_immuno %in% rcdf_assays, tp, 2])),
          xbreaks = seq(min(assay_lim[assay_immuno %in% rcdf_assays, tp, 1]), 
                        max(assay_lim[assay_immuno %in% rcdf_assays, tp, 2]), 
                        2),
          plot_title = paste0(
            switch(tp, Delta29overB = "Day 29", Delta57overB = "Day 57"),
            " Fold-rise over Baseline Ab Markers"
          ),
          filename = paste0(
            save.results.to, "/Marker_Rcdf_", c("nAb", "bAb")[bAb + 1], "_", tp,
            "_trt_vaccine_bstatus_", c("Neg", "Pos")[bstatus], "_", study_name, ".png"
          )
        )
      }
      
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
print("Boxplots 1:")

tps <- c("Day29", "Day57", "Delta29overB", "Delta57overB")
for (bstatus in 1:2) {
  for (tp in tps[tps %in% times]) {

    covid_corr_boxplot_facets(
      plot_dat = subset(
        dat.long.twophase.sample,
        Bserostatus == bstatus.labels[bstatus]
      ),
      x = "Trt",
      y = tp,
      color = "Trt",
      facet_by = "assay",
      ylim = assay_lim[, tp, ],
      plot_LLOD = tp %in% c("B", "Day29", "Day57"),
      LLOD = log10(llods[assay_immuno]),
      POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
      LLOQ = log10(lloqs[assay_immuno]),
      ULOQ = log10(uloqs[assay_immuno]),
      arrange_ncol = 3,
      arrange_nrow = ceiling(length(assay_immuno) / 3),
      legend = c("Placebo", "Vaccine"),
      axis_titles_y = labels.axis[tp, ] %>% unlist(),
      panel_titles = labels.title2[tp, ] %>% unlist(),
      filename = paste0(
        save.results.to, "/boxplots_", tp, "_x_trt_", bstatus.labels.2[bstatus],
        "_", study_name, ".png"
      )
    )
  }
}

#-----------------------------------------------
# - Box plots of the assay readouts versus baseline sero-status, stratified by
#   treatment groups
# - Make seperate plots for Placebo and Vaccine arms
#-----------------------------------------------
tps <- c("Day29", "Day57", "Delta29overB", "Delta57overB")
for (trt in 1:2) {
  for (tp in tps[tps %in% times]) {
    
    covid_corr_boxplot_facets(
      plot_dat = subset(dat.long.twophase.sample, as.numeric(Trt) == trt),
      x = "Bserostatus",
      y = tp,
      color = "Bserostatus",
      facet_by = "assay",
      ylim = assay_lim[, tp,],
      plot_LLOD = tp %in% c("B", "Day29", "Day57"),
      LLOD = log10(llods[assay_immuno]),
      POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
      LLOQ = log10(lloqs[assay_immuno]),
      ULOQ = log10(uloqs[assay_immuno]),
      arrange_ncol = 3,
      arrange_nrow = ceiling(length(assay_immuno) / 3),
      legend = c("Baseline Negative", "Baseline Positive"),
      axis_titles_y = labels.axis[tp, ] %>% unlist(),
      panel_titles = labels.title2[tp, ] %>% unlist(),
      filename = paste0(
        save.results.to, "/boxplots_", tp,
        "_x_bstatus_", c("placebo_arm_", "vaccine_arm_")[trt],
        study_name, ".png"
      )
    )
  }
}


#-----------------------------------------------
# Spaghetti PLOTS
#-----------------------------------------------
# - Spaghetti plots of antibody marker change over time
#-----------------------------------------------
print("Spaghetti plots:")
## in each baseline serostatus group, randomly select 10 placebo recipients and 20 vaccine recipients
set.seed(12345)
var_names <- expand.grid(times = intersect(c("B", "Day29", "Day57"), times),
                         assays = assay_immuno) %>%
  mutate(var_names = paste0(times, assay_immuno)) %>%
  .[, "var_names"]

spaghetti_ptid <- dat.twophase.sample[, c("Ptid", "Bserostatus", "Trt", var_names)] %>%
  filter(., complete.cases(.)) %>%
  transmute(BT = paste0(as.character(Bserostatus), as.character(Trt)),
            Ptid = Ptid) %>%
  split(., .$BT) %>%
  lapply(function(xx) {
    if (xx$BT[1] %in% c("10", "00")) {
      sample(xx$Ptid, 10)  ## sample 10 placebo recipients
    } else {
      sample(xx$Ptid, 20)  ## sample 20 vaccine recipients
    }
  }) %>% unlist %>% as.character

spaghetti_dat <- dat.long.twophase.sample[, c("Ptid", "Bserostatus", "Trt", "assay",
                                              intersect(c("B", "Day29", "Day57"), times))] %>%
  filter(Ptid %in% spaghetti_ptid) %>%
  pivot_longer(cols = intersect(c("B", "Day29", "Day57"), times),
               names_to = "time") %>%
  mutate(assay = factor(assay, levels = assay_immuno, labels = assay_immuno),
         time_label = factor(time, levels = intersect(c("B", "Day29", "Day57"), times),
                             labels = c("D1", "D29", "D57")[c("B", "Day29", "Day57") %in% times])) %>%
  as.data.frame

for (bstatus in 1:2) {
  subdat <- subset(spaghetti_dat, Bserostatus == bstatus.labels[bstatus])
  covid_corr_spaghetti_facets(plot_dat = subdat,
                              x = "time_label",
                              y = "value",
                              id = "Ptid",
                              color = "Trt",
                              facet_by = "assay",
                              ylim = assay_lim[, intersect(c("Day29", "Day57"), times)[1],],
                              panel_titles = labels.assays.short,
                              plot_title = paste0(
                                "Baseline ",
                                c("Negative", "Positive")[bstatus],
                                " PP Placebo + Vaccine group"
                              ),
                              arrange_ncol = 3,
                              arrange_nrow = ceiling(length(assay_immuno) / 3),
                              filename = paste0(
                                save.results.to, "/spaghetti_plot_",
                                bstatus.labels.2[bstatus], "_",
                                study_name, ".png"
                              ))
}



#### Scatter plot, assay vs. age in years, (Day 1) Day 29 Day 57
print("Scatter plots:")
tps <- c("B", "Day29", "Day57")
for (tp in tps[tps %in% times]) {
  for (bstatus in 1:2) {
    for (trt in 1:2) {
      subdat <- dat.long.twophase.sample %>%
        filter(Bserostatus == bstatus.labels[bstatus], Trt == trt.labels[trt])
      
      ## setting the range of the axes
      xrange <- range(dat.long.twophase.sample$Age)
      
      
      scatter_plot_list <- vector("list", length = length(assay_immuno))
      
      for (aa in 1:length(assay_immuno)) {
        scatter_plot_list[[aa]] <- ggplot(data = subset(subdat, assay == assay_immuno[aa]),
                                          mapping = aes_string("Age", tp)) +
          geom_point() +
          xlab("Age") +
          ylab(labels.axis[tp, aa]) +
          ggtitle(labels.title[tp, aa]) +
          stat_smooth(method = "loess", color = "red", se = TRUE, lwd = 1) +
          scale_x_continuous(limits = xrange) +
          scale_y_continuous(
            labels = label_math(10^.x), limits = assay_lim[aa, tp,],
            breaks = seq(assay_lim[aa, tp, 1], assay_lim[aa, tp, 2], by = 2)
          ) +
          theme_pubr() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 10),
            panel.border = element_rect(fill = NA),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_line(),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 10),
            legend.title = element_blank()
          )
      }
      
      output_plot <- ggarrange(
        plotlist = scatter_plot_list, ncol = 3,
        nrow = ceiling(length(assay_immuno) / 3), legend = "none", align = "h"
      )
      
      ggsave(
        filename = paste0(
          save.results.to, "/scatter_", tp, "_vs_age_",
          "trt_", trt.labels[trt], "_", bstatus.labels.2[bstatus],
          "_", study_name, ".png"
        ), 
        plot = output_plot, 
        width = 9,
        height = 0.5 + 3 * ceiling(length(assay_immuno) / 3), 
        units = "in"
      )
      
    }
  }
}


