source("descriptive_graphics_data_preprocess.R")
source("ggally_cor_resample.R")
source("covid_corr_plot_functions.R")

## =================================================================================================
## boxplots and weighted rcdf plots of assay readouts at different time points versus 
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
## =================================================================================================


for (tt in 1:5) {
  for (trt in 1:2) {
    for (bstatus in 1:2) {
      subdat <- subset(dat.long.twophase.sample,  
                       Bserostatus == bstatus.labels[bstatus] & 
                         Trt == trt.labels[trt])
      
      ##  (1) age >= 65 or age < 65 
      covid_corr_boxplot_by_assay(plot_dat = subdat,
                                  x = "age_geq_65_label", 
                                  y = times[tt],
                                  assays = assays, 
                                  plot_LLOQ = (tt <= 3), 
                                  LLOQ = LLOQ,
                                  axis_titles_y = labels.axis[tt, ],
                                  panel_titles = labels.title2[tt, ],
                                  filename = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],"_trt_", 
                                                    trt.labels[trt],"_by_age_group_", study.name,".png"))
      
      covid_corr_rcdf_by_assay(plot_dat = subdat,
                               time = times[tt],
                               assays = assays,
                               by = "age_geq_65_label",
                               weight = "wt",
                               panel_titles = labels.title2[tt, ],
                               axis_titles = labels.axis[tt, ],
                               filename = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                                                 "_trt_", c("placebo_", "vaccine_")[trt], bstatus.labels.2[bstatus], 
                                                 "_by_age_group_", study.name, ".png"))
      
      ##  (2) at risk / not at risk
      covid_corr_boxplot_by_assay(plot_dat = subdat,
                                  x = "highrisk_label", 
                                  y = times[tt],
                                  assays = assays, 
                                  plot_LLOQ = (tt <= 3), 
                                  LLOQ = LLOQ,
                                  axis_titles_y = labels.axis[tt, ],
                                  panel_titles = labels.title2[tt, ],
                                  filename = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],"_trt_", 
                                                    trt.labels[trt],"_by_risk_group_", study.name, ".png"))
      
      covid_corr_rcdf_by_assay(plot_dat = subdat,
                               time = times[tt],
                               assays = assays,
                               by = "highrisk_label",
                               weight = "wt",
                               panel_titles = labels.title2[tt, ],
                               axis_titles = labels.axis[tt, ],
                               filename = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                                                 "_trt_", c("placebo_", "vaccine_")[trt], bstatus.labels.2[bstatus], 
                                                 "_by_risk_group_", study.name, ".png"))
      
      ##  (3) age * high risk
      covid_corr_boxplot_by_assay(plot_dat = subdat,
                                  x = "age_risk_label", 
                                  y = times[tt],
                                  assays = assays,
                                  plot_LLOQ = (tt <= 3), 
                                  LLOQ = LLOQ,
                                  axis_titles_y = labels.axis[tt, ],
                                  panel_titles = labels.title2[tt, ],
                                  height = 7,
                                  filename = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],"_trt_", 
                                                    trt.labels[trt],"_by_age_x_risk_group_", study.name, ".png"))
      
      covid_corr_rcdf_by_assay(plot_dat = subdat,
                               time = times[tt],
                               assays = assays,
                               by = "age_risk_label",
                               weight = "wt",
                               panel_titles = labels.title2[tt, ],
                               axis_titles = labels.axis[tt, ],
                               height = 7,
                               filename = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                                                 "_trt_", c("placebo_", "vaccine_")[trt], bstatus.labels.2[bstatus], 
                                                 "_by_age_risk_group_", study.name, ".png"))
      
      ##  (4) sex at birth
      covid_corr_boxplot_by_assay(plot_dat = subdat,
                                  x = "sex_label", 
                                  y = times[tt],
                                  assays = assays, 
                                  plot_LLOQ = (tt <= 3), 
                                  LLOQ = LLOQ,
                                  axis_titles_y = labels.axis[tt, ],
                                  panel_titles = labels.title2[tt, ],
                                  filename = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],"_trt_", 
                                                    trt.labels[trt],"_by_sex_", study.name, ".png"))
      
      
      covid_corr_rcdf_by_assay(plot_dat = subdat,
                               time = times[tt],
                               assays = assays,
                               by = "sex_label",
                               weight = "wt",
                               panel_titles = labels.title2[tt, ],
                               axis_titles = labels.axis[tt, ],
                               filename = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                                                 "_trt_", c("placebo_", "vaccine_")[trt], bstatus.labels.2[bstatus], 
                                                 "_by_sex_group_", study.name, ".png"))
      
      ##  (5) age * sex at birth
      covid_corr_boxplot_by_assay(plot_dat = subdat,
                                  x = "age_sex_label", 
                                  y = times[tt],
                                  assays = assays,
                                  plot_LLOQ = (tt <= 3), 
                                  LLOQ = LLOQ,
                                  axis_titles_y = labels.axis[tt, ],
                                  panel_titles = labels.title2[tt, ],
                                  height = 7,
                                  filename = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],
                                                    "_trt_", trt.labels[trt],"_by_age_x_sex_group_", study.name, ".png"))
      
      covid_corr_rcdf_by_assay(plot_dat = subdat,
                               time = times[tt],
                               assays = assays,
                               by = "age_sex_label",
                               weight = "wt",
                               panel_titles = labels.title2[tt, ],
                               axis_titles = labels.axis[tt, ],
                               height = 7,
                               filename = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                                                 "_trt_", c("placebo_", "vaccine_")[trt], bstatus.labels.2[bstatus], 
                                                 "_by_age_sex_group_", study.name, ".png"))
      
      
      ##  (6) ethnicity
      covid_corr_boxplot_by_assay(plot_dat = subdat,
                                  x = "ethnicity", 
                                  y = times[tt],
                                  assays = assays,
                                  plot_LLOQ = (tt <= 3), 
                                  LLOQ = LLOQ,
                                  axis_titles_y = labels.axis[tt, ],
                                  panel_titles = labels.title2[tt, ],
                                  height = 7,
                                  filename = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],
                                                    "_trt_", trt.labels[trt],"_by_ethnicity_", study.name, ".png"))
      
      covid_corr_rcdf_by_assay(plot_dat = subdat,
                               time = times[tt],
                               assays = assays,
                               by = "ethnicity",
                               weight = "wt",
                               panel_titles = labels.title2[tt, ],
                               axis_titles = labels.axis[tt, ],
                               height = 7,
                               filename = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                                                 "_trt_", c("placebo_", "vaccine_")[trt], bstatus.labels.2[bstatus], 
                                                 "_by_ethnicity_", study.name, ".png"))
      
      ##  (7) race
      covid_corr_boxplot_by_assay(plot_dat = subset(subdat,
                                                    !(race == "White" & WhiteNonHispanic == 0)),
                                  x = "race", 
                                  y = times[tt],
                                  assays = assays,
                                  plot_LLOQ = (tt <= 3), 
                                  LLOQ = LLOQ,
                                  axis_titles_y = labels.axis[tt, ],
                                  panel_titles = labels.title2[tt, ],
                                  height = 7.5,
                                  filename = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],
                                                    "_trt_", trt.labels[trt],"_by_race_", study.name, ".png"))
      
      covid_corr_rcdf_by_assay(plot_dat = subset(subdat,
                                                 !(race == "White" & WhiteNonHispanic == 0)),
                               time = times[tt],
                               assays = assays,
                               by = "race",
                               weight = "wt",
                               panel_titles = labels.title2[tt, ],
                               axis_titles = labels.axis[tt, ],
                               height = 7.5,
                               filename = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                                                 "_trt_", c("placebo_", "vaccine_")[trt], bstatus.labels.2[bstatus], 
                                                 "_by_race_", study.name, ".png"))
      
      ##  (8) minority status
      covid_corr_boxplot_by_assay(plot_dat = subdat,
                                  x = "minority_label", 
                                  y = times[tt],
                                  assays = assays, 
                                  plot_LLOQ = (tt <= 3), 
                                  LLOQ = LLOQ,
                                  axis_titles_y = labels.axis[tt, ],
                                  panel_titles = labels.title2[tt, ],
                                  filename = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],
                                                    "_trt_", trt.labels[trt],"_by_minority_group_", study.name,".png"))
      
      covid_corr_rcdf_by_assay(plot_dat = subdat,
                               time = times[tt],
                               assays = assays,
                               by = "minority_label",
                               weight = "wt",
                               panel_titles = labels.title2[tt, ],
                               axis_titles = labels.axis[tt, ],
                               filename = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                                                 "_trt_", c("placebo_", "vaccine_")[trt], bstatus.labels.2[bstatus], 
                                                 "_by_minority_group_", study.name, ".png"))
      
      
      ##  (9) age * minority status 
      covid_corr_boxplot_by_assay(plot_dat = subdat,
                                  x = "age_minority_label", 
                                  y = times[tt],
                                  assays = assays,
                                  plot_LLOQ = (tt <= 3), 
                                  LLOQ = LLOQ,
                                  axis_titles_y = labels.axis[tt, ],
                                  panel_titles = labels.title2[tt, ],
                                  height = 7,
                                  filename = paste0(save.results.to, "demographics/boxplots_", times[tt], "_", bstatus.labels.2[bstatus],
                                                    "_trt_", trt.labels[trt],"_by_age_x_minority_", study.name, ".png"))
      
      covid_corr_rcdf_by_assay(plot_dat = subdat,
                               time = times[tt],
                               assays = assays,
                               by = "age_minority_label",
                               weight = "wt",
                               panel_titles = labels.title2[tt, ],
                               axis_titles = labels.axis[tt, ],
                               height = 7,
                               filename = paste0(save.results.to, "demographics/Marker_Rcdf_", times[tt], 
                                                 "_trt_", c("placebo_", "vaccine_")[trt], bstatus.labels.2[bstatus], 
                                                 "_by_age_minority_group_", study.name, ".png"))
      
    }
  }
}


