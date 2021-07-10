#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(dplyr)
library(stringr)
library(GGally)
library(ggpubr)
library(SWIM)
library(ggplot2)
library(scales)

set.seed(12345)
# load helper scripts and parameters
source(here("code", "ggally_cor_resample.R"))
source(here("code", "covid_corr_plot_functions.R"))
source(here("code", "params.R"))



# load cleaned data
dat.long.twophase.sample <- readRDS(here(
  "data_clean",
  "long_twophase_data.rds"
))
dat.twophase.sample <- readRDS(here("data_clean", "twophase_data.rds"))

assay_lim <- readRDS(here("data_clean", "assay_lim.rds"))
## ============================================================================
## boxplots and weighted rcdf plots of assay readouts at time points versus
##  (1) age >= 65 or age < 65
##  (2) at risk / not at risk
##  (3) age * high risk
##  (4) sex at birth
##  (5) age * sex at birth
##  (6) ethnicity
##  (7) race
##  (8) minority status
##  (9) age * minority status
## plot for each treatment group by baseline status
## ============================================================================

tps <- c("Day29", "Day57", "Delta29overB", "Delta57overB")
for (tp in tps[tps %in% times]) {
  for (trt in 1:2) {
    # Don't produce figures for placebo baseline negative to improve build time
    if(trt==1) {bstatus.range <- 2} else {bstatus.range <- 1:2}
    for (bstatus in bstatus.range) {
      subdat <- subset(
        dat.long.twophase.sample,
        Bserostatus == bstatus.labels[bstatus] &
          Trt == trt.labels[trt]
      )

      ##  (1) age >= 65 or age < 65
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "age_geq_65_label",
        y = tp,
        facet_by = "assay",
        plot_LLOD = tp %in% c("B", "Day29", "Day57"),
        LLOD = log10(llods[assay_immuno]),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOQ = log10(lloqs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        ylim = assay_lim[assay_immuno, tp, ],
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_ncol = 3,
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_age_group_",
          study_name, ".png"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[assay_immuno, tp, ],
        color = "age_geq_65_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_age_group_", study_name,
          ".png"
        )
      )

      ##  (2) at risk / not at risk
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "highrisk_label",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[assay_immuno, tp, ],
        plot_LLOD = tp %in% c("B", "Day29", "Day57"),
        LLOD = log10(llods[assay_immuno]),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOQ = log10(lloqs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_risk_group_",
          study_name, ".png"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[assay_immuno, tp, ],
        color = "highrisk_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_risk_group_", study_name,
          ".png"
        )
      )

      ##  (3) age * high risk
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "age_risk_label",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[assay_immuno, tp, ],
        plot_LLOD = tp %in% c("B", "Day29", "Day57"),
        LLOD = log10(llods[assay_immuno]),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOQ = log10(lloqs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_age_x_risk_group_",
          study_name, ".png"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[assay_immuno, tp, ],
        color = "age_risk_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_age_risk_group_",
          study_name, ".png"
        )
      )

      ##  (4) sex at birth
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "sex_label",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[assay_immuno, tp, ],
        plot_LLOD = tp %in% c("B", "Day29", "Day57"),
        LLOD = log10(llods[assay_immuno]),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOQ = log10(lloqs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_sex_", study_name,
          ".png"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[assay_immuno, tp, ],
        color = "sex_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_sex_group_", study_name,
          ".png"
        )
      )

      ##  (5) age * sex at birth
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "age_sex_label",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[assay_immuno, tp, ],
        plot_LLOD = tp %in% c("B", "Day29", "Day57"),
        LLOD = log10(llods[assay_immuno]),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOQ = log10(lloqs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_age_x_sex_group_",
          study_name, ".png"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[assay_immuno, tp, ],
        color = "age_sex_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_age_sex_group_",
          study_name, ".png"
        )
      )

      # (6) ethnicity
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "ethnicity_label",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[assay_immuno, tp, ],
        plot_LLOD = tp %in% c("B", "Day29", "Day57"),
        LLOD = log10(llods[assay_immuno]),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOQ = log10(lloqs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_ethnicity_", study_name,
          ".png"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[assay_immuno, tp, ],
        color = "ethnicity_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_ethnicity_", study_name,
          ".png"
        )
      )

      # (7) race
      covid_corr_boxplot_facets(
        plot_dat =
          subset(subdat, !(race == "White" &
            WhiteNonHispanic == 0)),
        x = "race",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[assay_immuno, tp, ],
        plot_LLOD = tp %in% c("B", "Day29", "Day57"),
        LLOD = log10(llods[assay_immuno]),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOQ = log10(lloqs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_race_", study_name,
          ".png"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat =
          subset(subdat, !(race == "White" &
            WhiteNonHispanic == 0)),
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[assay_immuno, tp, ],
        color = "race",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_race_", study_name,
          ".png"
        )
      )

      ##  (8) minority status
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "minority_label",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[assay_immuno, tp, ],
        plot_LLOD = tp %in% c("B", "Day29", "Day57"),
        LLOD = log10(llods[assay_immuno]),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOQ = log10(lloqs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_minority_group_",
          study_name, ".png"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[assay_immuno, tp, ],
        color = "minority_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_minority_group_",
          study_name, ".png"
        )
      )

      ##  (9) age * minority status
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "age_minority_label",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[assay_immuno, tp, ],
        plot_LLOD = tp %in% c("B", "Day29", "Day57"),
        LLOD = log10(llods[assay_immuno]),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOQ = log10(lloqs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_age_x_minority_",
          study_name, ".png"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[assay_immuno, tp, ],
        color = "age_minority_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_age_minority_group_",
          study_name, ".png"
        )
      )
      print(paste0("Done with loop of ", tp, ", trt=",
                   trt," and bstatus=", bstatus))
    }
  }
}
