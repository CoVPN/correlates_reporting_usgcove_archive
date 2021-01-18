source("descriptive_graphics_data_preprocess.R")
source("ggally_cor_resample.R") ## produce geom_statistics with resampling-based covariate-adjusted Spearman correlation
####################################
###          PAIR PLOTS          ### 
####################################


#================================================================================================================================================================================
# The correlation of each pair of Day 57 and delta antibody marker readouts are compared  within each of the baseline strata subgroups, treatment arm, and baseline serostratum.  
# Pairs plots/scatterplots and baseline strata-adjusted Spearman rank correlations are used. 
# Again, we made the correlation plots for each serostratum and all serostrata combined
#================================================================================================================================================================================

# for (tt in 1:2) {
#   for (bserostatus in 0:1) {
#     for (trt in c(0,1)) {    
#       for (bstratum in 1:max.stratum) {
#         dat.tmp <- subset(dat.mock, Bserostatus == bserostatus & Bstratum2 == bstratum & Trt == trt)[, c("Day57", "Delta")[tt] %.% assays]
#         main <- paste0("D1 Ab markers: baseline ", ifelse(bserostatus, "positive", "negative"), "\n", 
#                        c("placebo", "vaccine")[trt + 1], " arm, stratum ", bstratum)
#         png(filename = paste0(save.results.to, "stratum/pairs_", c("Day57", "Delta")[tt],"_Markers_",
#                               c("placebo", "vaccine")[trt + 1], "_arm_", bstatus.labels.2[bserostatus+1],"_stratum_",bstratum, ".png"),
#             width = 1520, height = 1540)
#         pairplots <- ggpairs(data = dat.tmp, title = main,
#                              columnLabels = labels.axis[1 + tt, ], upper = list(continuous = wrap(ggally_cor, stars = FALSE, size = 18))) + theme_bw() +
#           theme(plot.title = element_text(hjust = 0.5, size = 28),
#                 strip.text = element_text(size = 20, face = "bold"),
#                 axis.text = element_text(size = 28),
#                 panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
#         pairplots[1, 1] <- pairplots[1, 1] + scale_x_continuous(limits = c(-1, 6), labels = label_math(10 ^ .x)) + ylim(0, 1.2)
#         for(j in 2:pairplots$nrow) {
#           for (k in 1:(j - 1)) {
#             pairplots[j, k] <- pairplots[j, k] + 
#               stat_smooth(method = "loess", color = "red", se = FALSE) + 
#               scale_x_continuous(limits = c(-1, 6), labels = label_math(10 ^ .x)) +
#               scale_y_continuous(limits = c(-1, 6), labels = label_math(10 ^ .x))
#           }
#           pairplots[j, j] <- pairplots[j, j] + scale_x_continuous(limits = c(-1, 6), labels = label_math(10 ^ .x)) + ylim(0, 1.2)
#         }
#         
#         
#         print(pairplots)
#         dev.off()
#         
#       }
#     }
#   }
# }
# 
# 
# 
# 
# for (tt in 1:2) {
#   for (bserostatus in 0:1) {
#     for (bstratum in 1:max.stratum) {
#       dat.tmp <- subset(dat.mock, Bserostatus == bserostatus & Bstratum2 == bstratum & Trt == trt)[, c("Day57", "Delta")[tt] %.% assays]
#       main <- paste0("D1 Ab markers: baseline ", ifelse(bserostatus, "positive", "negative"), ", placebo + vaccine arm, stratum ", bstratum)
#       png(filename = paste0(save.results.to, "stratum/pairs_", c("Day57", "Delta")[tt],"_Markers_",bstatus.labels.2[bserostatus+1],"_stratum",bstratum, ".png"),
#           width = 1520, height = 1540)
#       pairplots <- ggpairs(data = dat.tmp, title = main,
#                            columnLabels = labels.axis[1 + tt, ], upper = list(continuous = wrap(ggally_cor, stars = FALSE, size = 16))) + theme_bw() +
#         theme(plot.title = element_text(hjust = 0.5, size = 28),
#               strip.text = element_text(size = 20, face = "bold"),
#               axis.text = element_text(size = 28),
#               panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
#       pairplots[1, 1] <- pairplots[1, 1] + scale_x_continuous(limits = c(-1, 6), labels = label_math(10 ^ .x)) + ylim(0, 1.2)
#       for(j in 2:pairplots$nrow) {
#         for (k in 1:(j - 1)) {
#           pairplots[j, k] <- pairplots[j, k] + 
#             stat_smooth(method = "loess", color = "red", se = FALSE) + 
#             scale_x_continuous(limits = c(-1, 6), labels = label_math(10 ^ .x)) +
#             scale_y_continuous(limits = c(-1, 6), labels = label_math(10 ^ .x))
#         }
#         pairplots[j, j] <- pairplots[j, j] + scale_x_continuous(limits = c(-1, 6), labels = label_math(10 ^ .x)) + ylim(0, 1.2)
#       }
#       
#       
#       print(pairplots)
#       dev.off()
#       
#     }
#     
#   }
# }




for (tt in 1:4) {
  for (bserostatus in 0:1) {
    for (trt in c(0,1)) {    
      subdat <- subset(dat.mock.twophase.sample, Bserostatus==bserostatus & Trt==trt)
      dat.tmp <- subdat[, times[tt + 1] %.% assays]
      strata <- subdat$Bstratum
      weight <- subdat$wt
        
      main <- paste0(c("D29", "D57", "D29 Fold-rise over D1", "D57 Fold-rise over D1")[tt], " Ab markers: baseline ", ifelse(bserostatus, "positive", "negative"), ", ", 
                     c("placebo", "vaccine")[trt + 1], " arm")
      
      rr <- range(dat.tmp, na.rm = TRUE) 
      
      if (rr[2] - rr[1] < 2) {
        rr <- floor(rr[1]):ceiling(rr[2])
      }
      
      breaks <- floor(rr[1]):ceiling(rr[2])
      
      if (rr[2] > ceiling(rr[1])) {
        breaks <- ceiling(rr[1]):floor(rr[2])
      } else {
        breaks <- floor(rr[1]):ceiling(rr[2]) ## breaks on the axis
      }
      
      if (max(breaks) - min(breaks) >= 6) {
        breaks <- breaks[breaks %% 2 == 0]
      }
      
      
      
      pairplots <- ggpairs(data = dat.tmp, title = main,
                           columnLabels = labels.axis[tt + 1, 1:4], 
                           upper = list(continuous = wrap(ggally_cor_resample, 
                                                          stars = FALSE, 
                                                          size = 16, 
                                                          strata = strata,
                                                          weight = weight))) + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size = 36),
              strip.text = element_text(size = 28, face = "bold"),
              axis.text = element_text(size = 28),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
      pairplots[1, 1] <- pairplots[1, 1] + scale_x_continuous(limits = rr, breaks = breaks) + ylim(0, 1.2)
      for(j in 2:pairplots$nrow) {
        for (k in 1:(j - 1)) {
          pairplots[j, k] <- pairplots[j, k] + 
            stat_smooth(method = "loess", color = "red", se = FALSE, lwd = 2) + 
            scale_x_continuous(limits = rr, breaks = breaks, labels = label_math(10 ^ .x)) +
            scale_y_continuous(limits = rr, breaks = breaks, labels = label_math(10 ^ .x))
        }
        pairplots[j, j] <- pairplots[j, j] + scale_x_continuous(limits = rr, breaks = breaks, labels = label_math(10 ^ .x)) + ylim(0, 1.2)
      }
      
      png(filename = paste0(save.results.to, "pairs_", times[tt + 1], 
                            "_Markers_", bstatus.labels.2[bserostatus+1], 
                            c("_placebo_arm", "_vaccine_arm")[trt + 1], "_",
                            study.name, ".png"),
          width = 1520, height = 1540)
      print(pairplots)
      dev.off()
      
    }
  }
}



for (tt in 1:4) {
  for (bserostatus in 0:1) {
    subdat <- subset(dat.mock.twophase.sample, Bserostatus==bserostatus)
    dat.tmp <- subdat[, times[tt + 1] %.% assays]
    strata <- subdat$Bstratum
    weight <- subdat$wt
    
    rr <- range(dat.tmp, na.rm = TRUE) 
    
    if (rr[2] - rr[1] < 2) {
      rr <- floor(rr[1]):ceiling(rr[2])
    }
    
    breaks <- floor(rr[1]):ceiling(rr[2])
    
    if (rr[2] > ceiling(rr[1])) {
      breaks <- ceiling(rr[1]):floor(rr[2])
    } else {
      breaks <- floor(rr[1]):ceiling(rr[2]) ## breaks on the axis
    }
    
    if (max(breaks) - min(breaks) >= 6) {
      breaks <- breaks[breaks %% 2 == 0]
    }
    
    
    main <- paste0(c("D29", "D57", "D29 Fold-rise over D1", "D57 Fold-rise over D1")[tt], " Ab markers: baseline ", ifelse(bserostatus, "positive", "negative"), 
                   ", placebo + vaccine arm")
    png(filename = paste0(save.results.to, "pairs_", c("Day29", "Day57", "Delta29overB", "Delta57overB")[tt], 
                          "_Markers_", bstatus.labels.2[bserostatus+1], "_", study.name, ".png"),
        width = 1520, height = 1540)
    pairplots <- ggpairs(data = dat.tmp, title = main,
                         columnLabels = labels.axis[tt + 1,], 
                         upper = list(continuous = wrap(ggally_cor_resample, 
                                                        stars = FALSE, 
                                                        strata = strata,
                                                        weight = weight,
                                                        size = 16))) + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 36),
            strip.text = element_text(size = 28, face = "bold"),
            axis.text = element_text(size = 28),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    pairplots[1, 1] <- pairplots[1, 1] + scale_x_continuous(limits = rr) + ylim(0, 1.2)
    for(j in 2:pairplots$nrow) {
      for (k in 1:(j - 1)) {
        pairplots[j, k] <- pairplots[j, k] + 
          stat_smooth(method = "loess", color = "red", se = FALSE, lwd = 2) +
          scale_x_continuous(limits = rr, breaks = breaks, labels = label_math(10 ^ .x)) +
          scale_y_continuous(limits = rr, breaks = breaks, labels = label_math(10 ^ .x))
      }
      pairplots[j, j] <- pairplots[j, j] + scale_x_continuous(limits = rr, breaks = breaks, labels = label_math(10 ^ .x)) + ylim(0, 1.2)
    }
    
    
    print(pairplots)
    dev.off()
    
  }
  
}


#===========================================================================================
# The correlation of each pair of Day 1 antibody marker readouts are compared within each of
# the baseline demographic subgroups and baseline serostratum, pooling over the two treatment
# arms. Pairs plots/scatterplots and Spearman rank correlations are used, where the plots and
# tests pooling over the baseline demographic subgroups use the case-deleted data set.


for (bserostatus in 0:1) {
  for (trt in 0:1) { 
    subdat <- subset(dat.mock.twophase.sample, Bserostatus == bserostatus & as.numeric(Trt) == trt)
    dat.tmp <- subdat[, "B"%.%assays]
    strata <- subdat$Bstratum
    weight <- subdat$wt
    
    rr <- range(dat.tmp, na.rm = TRUE) 
    
    if (rr[2] - rr[1] < 2) {
      rr <- floor(rr[1]):ceiling(rr[2])
    }
    
    breaks <- floor(rr[1]):ceiling(rr[2])
    
    if (rr[2] > ceiling(rr[1])) {
      breaks <- ceiling(rr[1]):floor(rr[2])
    } else {
      breaks <- floor(rr[1]):ceiling(rr[2]) ## breaks on the axis
    }
    
    if (max(breaks) - min(breaks) >= 6) {
      breaks <- breaks[breaks %% 2 == 0]
    }
    
    png(filename = paste0(save.results.to, "pairs_baselineMarkers_",bstatus.labels.2[bserostatus + 1],
                          "_", c("placebo", "vaccine")[trt + 1], "_arm_", study.name, ".png"),
        width = 1520, height = 1540)
    main <- paste0("D1 Ab markers: baseline ", ifelse(bserostatus, "positive", "negative"), ", ", 
                   c("placebo", "vaccine")[trt + 1], " arm")
    pairplots <- ggpairs(data = dat.tmp, columnLabels = labels.axis[1, ], 
                         upper = list(continuous = wrap(ggally_cor_resample, 
                                                        stars = FALSE, 
                                                        size = 16, 
                                                        method = "spearman",
                                                        strata = strata,
                                                        weight = weight))) + theme_bw() +
      ggtitle(main) +
      theme(plot.title = element_text(hjust = 0.5, size = 36),
            strip.text = element_text(size = 28, face = "bold"),
            axis.text = element_text(size = 28),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    pairplots[1, 1] <- pairplots[1, 1] + scale_x_continuous(limits = rr, breaks = breaks) + ylim(0, 1.2)
    for(j in 2:pairplots$nrow) {
      for (k in 1:(j - 1)) {
        pairplots[j, k] <- pairplots[j, k] + 
          stat_smooth(method = "loess", color = "red", se = FALSE, lwd = 2) +
          scale_x_continuous(labels = label_math(10^.x), limits = rr, breaks = breaks) +
          scale_y_continuous(limits = rr, labels = label_math(10^.x), breaks = breaks)
      }
      pairplots[j, j] <- pairplots[j, j] + scale_x_continuous(limits = rr, labels = label_math(10^.x), breaks = breaks) + ylim(0, 1.2)
    }
    
    
    print(pairplots)
    dev.off()
  }
}




for (bserostatus in 0:1) {
  subdat <- subset(dat.mock.twophase.sample, Bserostatus == bserostatus & as.numeric(Trt) == trt)
  dat.tmp <- subdat[, "B"%.%assays]
  strata <- subdat$Bstratum
  weight <- subdat$wt
  
  rr <- range(dat.tmp, na.rm = TRUE) 
  
  if (rr[2] - rr[1] < 2) {
    rr <- floor(rr[1]):ceiling(rr[2])
  }
  
  breaks <- floor(rr[1]):ceiling(rr[2])
  
  if (rr[2] > ceiling(rr[1])) {
    breaks <- ceiling(rr[1]):floor(rr[2])
  } else {
    breaks <- floor(rr[1]):ceiling(rr[2]) ## breaks on the axis
  }
  
  if (max(breaks) - min(breaks) >= 6) {
    breaks <- breaks[breaks %% 2 == 0]
  }
  
  
  png(filename = paste0(save.results.to, "pairs_baselineMarkers_",bstatus.labels.2[bserostatus + 1],
                        "_", study.name, ".png"),
      width = 1520, height = 1540)
  main <- paste0("D1 Ab markers: baseline ", ifelse(bserostatus, "positive", "negative"), ", vaccine + placebo arm")
  pairplots <- ggpairs(data = dat.tmp, columnLabels = labels.axis[1, ], 
                       upper = list(continuous = wrap(ggally_cor_resample, 
                                                      stars = FALSE, 
                                                      size = 16, 
                                                      method = "spearman",
                                                      strata = strata,
                                                      weight = weight))) + 
    theme_bw() +
    ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5, size = 36),
          strip.text = element_text(size = 28, face = "bold"),
          axis.text = element_text(size = 28),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  pairplots[1, 1] <- pairplots[1, 1] + scale_x_continuous(limits = rr, breaks = breaks) + ylim(0, 1.2)
  for(j in 2:pairplots$nrow) {
    for (k in 1:(j - 1)) {
      pairplots[j, k] <- pairplots[j, k] + 
        stat_smooth(method = "loess", color = "red", se = FALSE, lwd = 2) +
        scale_x_continuous(labels = label_math(10^.x), rr, breaks = breaks) +
        scale_y_continuous(limits = rr, labels = label_math(10^.x), breaks = breaks)
    }
    pairplots[j, j] <- pairplots[j, j] + scale_x_continuous(limits = rr, labels = label_math(10^.x), breaks = breaks) + ylim(0, 1.2)
  }
  
  
  print(pairplots)
  dev.off()
  
}

####################################
###          RCDF PLOTS          ### 
####################################

## original plot size: 960 * 960

#=========================================================================================================================
# Reverse empirical cdf (rcdf) plots for the Baseline/Day 57/Baseline-subtracted Day 57assay readouts, stratified by treatment group and baseline serostatus
# We made four ggplot objects, each for one assay, and combine them with ggarrange
#=========================================================================================================================

for (tt in 1:5){
  rcdf_list <- vector("list", 4)
  for (aa in 1:4) {
    rcdf_list[[aa]] <- ggplot(subset(dat.mock.long.twophase.sample, assay == assays[aa]), aes_string(x = times[tt], colour = "Trt:Bserostatus", weight = "wt")) +
      geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr(legend = "none") + 
      ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) +
      scale_x_continuous(labels = label_math(10^.x), limits = c(-2, 10), breaks = seq(-2, 10, 2)) +
      scale_color_manual(values = hvtn_col,
                         labels = c("Placebo, Baseline Neg",
                                    "Placebo, Baseline Pos",
                                    "Vaccine, Baseline Neg",
                                    "Vaccine, Baseline Pos")) +
      ggtitle(labels.title2[tt, aa]) +
      theme(plot.title = element_text(hjust = 0.5, size = 24),
            legend.title = element_blank(),
            legend.text = element_text(size = 18, face = "bold"),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_line(),
            axis.title = element_text(size = 24),
            axis.text = element_text(size = 28))
  }
  png(file = paste0(save.results.to, "Marker_RCDF_", times[tt], "_trt_both_bstatus_both_", study.name, ".png"), 
      height = 1060, width = 1060)
  print(ggarrange(plotlist = rcdf_list, ncol = 2, nrow = 2, 
                  common.legend = TRUE, legend = "bottom",
                  align = "h"))
  dev.off()
}


#=========================================================================================================================
# Reverse empirical cdf (rcdf) plots for the Baseline/Day 57/Baseline-subtracted Day 57assay readouts, 
# stratified by treatment group and event status, in baseline negative or positive subjects
# We made four ggplot objects, each for one assay, and combine them with ggarrange
#=========================================================================================================================

for (bb in 1:2) {
  for (tt in c(1, 2, 4)){
    rcdf_list <- vector("list", 4)
    for (aa in 1:4) {
      rcdf_list[[aa]] <- ggplot(subset(dat.mock.long.twophase.sample, assay == assays[aa] & Bserostatus == bstatus.labels[bb]), 
                                aes_string(x = times[tt], colour = "Trt:EventLabelD29", linetype = "EventLabelD29", weight = "wt")) +
        geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr(legend = "none") +
        ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) +
        scale_x_continuous(labels = label_math(10^.x), limits = c(-2, 10), breaks = seq(-2, 10, 2)) +
        scale_color_manual(values = hvtn_col) +
        guides(linetype = "none") +
        ggtitle(labels.title2[tt, aa]) +
        theme(plot.title = element_text(hjust = 0.5, size = 24),
              legend.title = element_blank(),
              legend.text = element_text(size = 18, face = "bold"),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 24),
              axis.text = element_text(size = 28))
    }
    png(file = paste0(save.results.to, "Marker_RCDF_", times[tt], "_trt_by_caseD29_bstatus_", bstatus.labels.2[bb],
                      "_", study.name,".png"), 
        height = 1060, width = 1060)
    print(ggarrange(plotlist = rcdf_list, ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom",
                    align = "h"))
    dev.off()
  }
  
}


for (bb in 1:2) {
  for (tt in 1:5){
    rcdf_list <- vector("list", 4)
    for (aa in 1:4) {
      rcdf_list[[aa]] <- ggplot(subset(dat.mock.long.twophase.sample, assay == assays[aa] & Bserostatus == bstatus.labels[bb]), 
                                aes_string(x = times[tt], colour = "Trt:EventLabelD57", linetype = "EventLabelD57", weight = "wt")) +
        geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr(legend = "none") +
        ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) +
        scale_x_continuous(labels = label_math(10^.x), limits = c(-2, 10), breaks = seq(-2, 10, 2)) +
        scale_color_manual(values = hvtn_col) +
        guides(linetype = "none") +
        ggtitle(labels.title2[tt, aa]) +
        theme(plot.title = element_text(hjust = 0.5, size = 24),
              legend.title = element_blank(),
              legend.text = element_text(size = 18, face = "bold"),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 24),
              axis.text = element_text(size = 28))
    }
    png(file = paste0(save.results.to, "Marker_RCDF_", times[tt], "_trt_by_caseD57_bstatus_", bstatus.labels.2[bb],
                      "_", study.name,".png"), 
        height = 1060, width = 1060)
    print(ggarrange(plotlist = rcdf_list, ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom",
                    align = "h"))
    dev.off()
  }
  
}

#=========================================================================================================================
# Reverse empirical cdf (rcdf) plots for the Baseline/Day 57/Baseline-subtracted Day 57assay readouts, stratified by treatment group and baseline serostatus
# We made four ggplot objects, each for one assay, and combine them with ggarrange
#=========================================================================================================================

for (tt in 1:5){
  rcdf_list <- vector("list", 4)
  for (aa in 1:4) {
    rcdf_list[[aa]] <- ggplot(subset(dat.mock.long.twophase.sample, assay == assays[aa]), 
                              aes_string(x = times[tt], colour = "Trt:Bserostatus", weight = "wt")) +
      geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr(legend = "none") +
      ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) +
      scale_x_continuous(labels = label_math(10 ^.x), limits = c(-2, 10), breaks = seq(-2, 10, 2)) +
      scale_color_manual(values = hvtn_col,
                         labels = c("Placebo, Baseline Neg",
                                    "Placebo, Baseline Pos",
                                    "Vaccine, Baseline Neg",
                                    "Vaccine, Baseline Pos")) +
      ggtitle(labels.title2[tt, aa]) +
      theme(plot.title = element_text(hjust = 0.5, size = 24),
            legend.title = element_blank(),
            legend.text = element_text(size = 18, face = "bold"),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_line(),
            axis.title = element_text(size = 24),
            axis.text = element_text(size = 28))
  }
  png(file = paste0(save.results.to, "Marker_RCDF_", times[tt], "_trt_both_bstatus_both",
                    "_", study.name, ".png"), 
      height = 1060, width = 1060)
  print(ggarrange(plotlist = rcdf_list, ncol = 2, nrow = 2, 
                  common.legend = TRUE, legend = "bottom",
                  align = "h"))
  dev.off()
}
##==========================================================================
## RCDF plot of four day 29 assay readouts in one plot, with the line-types 
## distinguishing the baseline serostatus
##==========================================================================

png(file = paste0(save.results.to, "Marker_Rcdf_day29_trt_vaccine_bstatus_both_", study.name, ".png"), 
    height = 640, width = 1020)
ggplot(subset(dat.mock.long.v2, Trt == "Vaccine"), aes(x = Day29, colour = assay, lty = Bserostatus, weight = "wt")) +
  geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
  scale_x_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
  scale_color_manual(values = hvtn_col) +
  ylab("Reverse ECDF") + xlab(bquote("D29 Binding Ab (IU/ml) / Pseudovirus nAb ID50 or ID80")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        panel.grid.minor.y = element_line(),
        panel.grid.major.y = element_line(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 28))
dev.off()


##==========================================================================
## RCDF plot of four day 57 assay readouts in one plot, with the line-types 
## distinguishing the baseline serostatus
##==========================================================================

png(file = paste0(save.results.to, "Marker_Rcdf_day57_trt_vaccine_bstatus_both_",
                  study.name, ".png"), height = 640, width = 1020)
ggplot(subset(dat.mock.long.v2, Trt == "Vaccine"), aes(x = Day57, colour = assay, lty = Bserostatus, weight = "wt")) +
  geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
  scale_x_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
  scale_color_manual(values = hvtn_col) +
  ylab("Reverse ECDF") + xlab(bquote("D57 Binding Ab (IU/ml) / Pseudovirus nAb ID50 or ID80")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        panel.grid.minor.y = element_line(),
        panel.grid.major.y = element_line(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 28))
dev.off()

##==========================================================================
## RCDF plot of four day 29 / day 57 fold-rise over baseline assay readouts 
## in one plot, with the line-types distinguishing the baseline serostatus
##==========================================================================


png(file = paste0(save.results.to, "Marker_Rcdf_Delta29overB_trt_vaccine_bstatus_both_", 
                  study.name,".png"), height = 640, width = 1020)
ggplot(subset(dat.mock.long.v2, Trt == "Vaccine"), 
       aes(x = Delta29overB, colour = assay, lty = Bserostatus, weight = "wt")) +
  geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
  scale_x_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
  scale_color_manual(values = hvtn_col) +
  ylab("Reverse ECDF") + xlab(bquote("Day 29 Fold-rise over Day 1 Binding Ab (IU/ml) / Pseudo nAb ID50 or ID80")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        panel.grid.minor.y = element_line(),
        panel.grid.major.y = element_line(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 28)) 
dev.off()

png(file = paste0(save.results.to, "Marker_Rcdf_Delta57overB_trt_vaccine_bstatus_both_", study.name, ".png"), 
    height = 640, width = 1020)
ggplot(subset(dat.mock.long.v2, Trt == "Vaccine"), aes(x = Delta57overB, colour = assay, lty = Bserostatus, weight = "wt")) +
  geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
  scale_x_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
  scale_color_manual(values = hvtn_col) +
  ylab("Reverse ECDF") + xlab(bquote("Day 57 Fold-rise over Day 1 Binding Ab (IU/ml) / Pseudo nAb ID50 or ID80")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        panel.grid.minor.y = element_line(),
        panel.grid.major.y = element_line(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 28)) 
dev.off()

##==========================================================================
## RCDF plot of four day 57 assay readouts in one plot for only the vaccine
## recipients with baseline negative sero viral status
##==========================================================================

png(file = paste0(save.results.to, "Marker_Rcdf_Day29_trt_vaccine_bstatus_Neg_", study.name,".png"), 
    height = 640, width = 1020)
ggplot(subset(dat.mock.long.v2, Trt == "Vaccine" & Bserostatus == "Baseline Neg"), aes(x = Day29, colour = assay, weight = "wt")) +
  geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
  scale_x_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
  scale_color_manual(values = hvtn_col) +
  ylab("Reverse ECDF") + xlab(bquote("D29 Binding Ab (IU/ml) / Pseudo nAb ID50 or ID80")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        panel.grid.minor.y = element_line(),
        panel.grid.major.y = element_line(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 28)) 
dev.off()



png(file = paste0(save.results.to, "Marker_Rcdf_Day57_trt_vaccine_bstatus_Neg_", study.name,".png"), height = 640, width = 1020)
ggplot(subset(dat.mock.long.v2, Trt == "Vaccine" & Bserostatus == "Baseline Neg"), aes(x = Day57, colour = assay, weight = "wt")) +
  geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
  scale_x_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
  scale_color_manual(values = hvtn_col) +
  ylab("Reverse ECDF") + xlab(bquote("D57 Binding Ab (IU/ml) / Pseudo nAb ID50 or ID80")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        panel.grid.minor.y = element_line(),
        panel.grid.major.y = element_line(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 28)) 
dev.off()




png(file = paste0(save.results.to, "Marker_Rcdf_Delta29_trt_vaccine_bstatus_Neg_", study.name, ".png"), 
    height = 640, width = 1020)
ggplot(subset(dat.mock.long.v2, Trt == "Vaccine" & Bserostatus == "Baseline Neg"), aes(x = Delta29overB, colour = assay, weight = "wt")) +
  geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
  scale_x_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
  scale_color_manual(values = hvtn_col) +
  ylab("Reverse ECDF") + xlab(bquote("Day 29 Fold-rise over Day 1 Binding Ab (IU/ml) / Pseudo nAb ID50 or ID80")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        panel.grid.minor.y = element_line(),
        panel.grid.major.y = element_line(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20)) 
dev.off()


png(file = paste0(save.results.to, "Marker_Rcdf_Delta57_trt_vaccine_bstatus_Neg_", study.name,".png"), 
    height = 640, width = 1020)
ggplot(subset(dat.mock.long.v2, Trt == "Vaccine" & Bserostatus == "Baseline Neg"), aes(x = Delta57overB, colour = assay, weight = "wt")) +
  geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
  scale_x_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
  scale_color_manual(values = hvtn_col) +
  ylab("Reverse ECDF") + xlab(bquote("Day 57 Fold-rise over Day 1 Binding Ab (IU/ml) / Pseudo nAb ID50 or ID80")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        panel.grid.minor.y = element_line(),
        panel.grid.major.y = element_line(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20)) 
dev.off()

##======================================================================================
## for the baseline-negative vaccine recipients (main group of interest), make cdf plots
## for each day 57 assay readouts. the color of the curves represent different subgroups
## defined by age, minority and if they belong in the high-risk category.
## We make four ggplot subjects, each for one assay, and combine them with ggarrange.
##======================================================================================
# 
# sub.dat.long <- subset(dat.mock.long, Trt == "Vaccine" & Bserostatus == "Baseline Neg")
# 
# for(tt in 2:3) {
#   rcdf_list <- vector("list", 4)
#   for (aa in 1:4) {
#     rcdf_list[[aa]] <- ggplot(subset(sub.dat.long, assay == assays[aa]), aes_string(x = times2[tt], colour = "Bstratum2", weight = "wt")) +
#       geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1.2) +  theme_pubr() +
#       ggtitle(labels.title[tt, aa]) +
#       scale_x_continuous(limits = c(-1, 6), labels = label_math(10^.x)) +
#       scale_color_manual(values = hvtn_col,
#                          labels = bstrata.labels) +
#       ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) + 
#       guides(colour = guide_legend(byrow = TRUE, nrow = 3)) +
#       theme(plot.title = element_text(hjust = 0.5, size = 20),
#             panel.grid.minor.y = element_line(),
#             panel.grid.major.y = element_line(),
#             axis.title = element_text(size = 20),
#             axis.text = element_text(size = 28),
#             legend.title = element_blank(),
#             legend.text = element_text(size = 20)) 
#   }
#   
#   png(file = paste0(save.results.to, "Marker_Rcdf_", times2[tt], "_trt_vaccine_bstatus_Neg_stratified.png"), height = 1240, width = 1120)
#   print(ggarrange(plotlist = rcdf_list , ncol = 2, nrow = 2, 
#                   common.legend = TRUE, legend = "bottom", 
#                   align = "h") + theme(plot.title = element_text(hjust = 0.5, size = 24)))
#   dev.off()
# }
########################################
###           SCATTER PLOTS          ### 
########################################




##========================================================================================
## Scatter plots for correlation between Day29/Day57/Delta29overB/Delta57overB and baseline, stratified by treatment,
## assay type and baseline serostatus.
## We created four ggplot objects, each for one assay. 
## Each ggplot object is the scatter plot of Day 57 assay readout or increase/fold-raise 
## from baseline versus the baseline assay readouts, stratified by the treatment group.
##========================================================================================


for (tt in 2:5) {
  for (bserostatus in 1:2) {
    for (trt in 1:2) {
      subdat <- subset(dat.mock.long.twophase.sample, as.numeric(Bserostatus) == bserostatus &
                         as.numeric(Trt) == trt)[, c("B", "Day29", "Day57", "Delta29overB", 
                                                     "Delta57overB", "assay", "Trt", "Bstratum", "wt")]
      scatterplot_list <- vector("list", 4)
      
      ## make the plot axis limits adaptive to the data range
      rr <- range(subdat[, c("B", times[tt])], na.rm = TRUE) 
      
      if (rr[2] - rr[1] < 2) {
        rr <- floor(rr[1]):ceiling(rr[2])
      }
      
      breaks <- floor(rr[1]):ceiling(rr[2])
      
      if (rr[2] > ceiling(rr[1])) {
        breaks <- ceiling(rr[1]):floor(rr[2])
      } else {
        breaks <- floor(rr[1]):ceiling(rr[2]) ## breaks on the axis
      }
      
      if (max(breaks) - min(breaks) >= 6) {
        breaks <- breaks[breaks %% 2 == 0]
      }
      
      
      for (aa in 1:4) {
        ## correlation
        
        ss <- subset(subdat, assay == assays[aa]) %>% filter(complete.cases(.))
  

        marker_corr <- round(spearman_resample(x = ss[, "B"], y = ss[, times[tt]], 
                                               strata = ss$Bstratum, weight = ss$wt, B = 200),
                             2)
        
        
        scatterplot_list[[aa]] <- ggscatter(x = "B", y = times[tt], data = subset(subdat, assay == assays[aa]), 
                                            title = paste0(assay.labels[aa]),
                                            xlab = labels.axis[1, aa], ylab = labels.axis[tt, aa], size = 2) +
          stat_smooth(method = "loess", color = "red", se = FALSE, lwd = 2) +
          scale_x_continuous(labels = label_math(10^.x), limits = rr, breaks = breaks) + 
          scale_y_continuous(labels = label_math(10^.x), limits = rr, breaks = breaks)+ 
          geom_text(x = 0.85 * rr[2] + 0.15 * rr[1], y = 0.93 * rr[2] + 0.07 * rr[1], label = paste0("Cor: ", marker_corr), size = 28 * 0.36) +
          theme(plot.title = element_text(hjust = 0.5, size = 20),
                panel.border = element_rect(fill = NA),
                panel.grid.minor.y = element_line(),
                panel.grid.major.y = element_line(),
                axis.title = element_text(size = 20),
                axis.text = element_text(size = 24),
                strip.text = element_text(size = 20, face = "bold"),
                legend.title = element_blank()) 
      }
      png(file = paste0(save.results.to, "scatterplots_", times[tt], "vB_", bstatus.labels.2[bserostatus],
                        c("_placebo_arm_", "_vaccine_arm_")[trt], study.name, ".png"), 
          height = 1080, width = 1040)
      print(ggarrange(plotlist = scatterplot_list, ncol = 2, nrow = 2,
                      legend = "none", common.legend = FALSE,
                      align = "h") + 
              theme(plot.title = element_text(hjust = 0.5, size = 24)))
      dev.off()
    }
  }
}









####################################
###           BOX PLOTS          ### 
####################################

## Box plots across treatment groups
## For the box plots, we made a ggplot object for every assay and use the ggarrange function to combine the resulted plots.

## =================================================================================================
## boxplots of Day 57 assay readouts versus treatment groups, stratified by baseline sero status,
## pooled over baseline strata
## =================================================================================================

## make a subset of data with 30 sample points for the jitter in each subgroup defined by Trt:Bserostatus
dat.mock.sample1 <- dat.mock.long.twophase.sample[, c("Trt", "Bserostatus", "B", "Day29", "Day57", "Delta29overB", "Delta57overB", "assay")] %>%
  filter(., complete.cases(.)) %>% 
  split(., list(.$Trt, .$Bserostatus, .$assay)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows



for (tt in 1:5){
  boxplot_list <- vector("list", 4)
  for (aa in 1:4) {
    boxplot_list[[aa]] <- ggplot(subset(dat.mock.long.twophase.sample, assay == assays[aa]), aes_string(x = "Trt", y = times[tt])) +
      geom_boxplot(aes(colour = Trt:Bserostatus), width = 0.6, lwd = 1.5) + 
      stat_boxplot(aes(colour = Trt:Bserostatus), geom = "errorbar", width = 0.45, lwd = 1.5) +
      guides(alpha = "none", fill = "none") +
      geom_jitter(data = subset(dat.mock.sample1, assay == assays[aa]), mapping = aes(colour = Trt:Bserostatus), width = 0.1, size = 5) +
      facet_wrap(vars(Bserostatus), scales = "free", ncol = 1) + 
      scale_y_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
      theme_pubr(legend = "none") + 
      ylab(labels.axis[tt, aa]) + xlab("") + 
      scale_fill_manual(values = hvtn_col) +
      scale_color_manual(values = hvtn_col, labels = c("Placebo, Baseline Neg",
                                                       "Placebo, Baseline Pos",
                                                       "Vaccine, Baseline Neg",
                                                       "Vaccine, Baseline Pos")) +
      ggtitle(labels.title2[aa]) +
      theme(plot.title = element_text(hjust = 0.5, size = 24),
            panel.border = element_rect(fill = NA),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_line(),
            axis.title = element_text(size = 24),
            axis.text = element_text(size = 28),
            strip.text = element_text(size = 20, face = "bold"),
            legend.title = element_blank(),
            legend.text = element_text(size = 20, face = "bold"))
    
    if (tt <= 3) {
      boxplot_list[[aa]] <- boxplot_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
        geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 8) 
    }
  }
  png(file = paste0(save.results.to, "boxplots_", times[tt], "_x_trt_", study.name,".png"), 
      height = 1840, width = 1040)
  print(ggarrange(plotlist = boxplot_list, ncol = 2, nrow = 2, 
                  common.legend = TRUE, legend = "bottom",
                  align = "h")) 
  dev.off()
}





#================================================================================
# The same plots but only for baseline negative subjects (main group of interest)
#================================================================================

for (bstatus in 1:2) {
  for (tt in 1:5){
    boxplot_list <- vector("list", 4)
    for (aa in 1:4) {
      boxplot_list[[aa]] <- ggplot(subset(dat.mock.long.twophase.sample, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus]), 
                                   aes_string(x = "Trt", y = times[tt])) +
        geom_boxplot(aes(colour = Trt), width = 0.6, lwd = 1.5) + 
        stat_boxplot(aes(colour = Trt), geom = "errorbar", width = 0.45, lwd = 1.5) +
        guides(alpha = "none", fill = "none") +
        geom_jitter(data = subset(dat.mock.sample1, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus]), mapping = aes(colour = Trt), width = 0.1, size = 5) +
        scale_y_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
        theme_pubr(legend = "none") + 
        ylab(labels.axis[tt, aa]) + xlab("") + 
        scale_fill_manual(values = hvtn_col) +
        scale_color_manual(values = hvtn_col, labels = c("Placebo", "Vaccine"))+
        ggtitle(labels.title2[tt, aa]) +
        theme(plot.title = element_text(hjust = 0.5, size = 24),
              panel.border = element_rect(fill = NA),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 28),
              strip.text = element_text(size = 20, face = "bold"),
              legend.title = element_blank(),
              legend.text = element_text(size = 24, face = "bold"))
      
      if (tt <= 3) {
        boxplot_list[[aa]] <- boxplot_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
          geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 8) 
      }
    }
    png(file = paste0(save.results.to, "boxplots_", times[tt], "_x_trt_", bstatus.labels.2[bstatus],
                      "_", study.name, ".png"), height = 1260, width = 1260)
    print(ggarrange(plotlist = boxplot_list, ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom",
                    align = "h")) 
    dev.off()
  }
}








## ==============================================================================================
## Box plots of the assay readouts versus baseline sero-status, stratified by treatment groups
## ==============================================================================================

for (tt in 1:5){
  boxplot_list <- vector("list", 4)
  for (aa in 1:4) {
    boxplot_list[[aa]] <- ggplot(subset(dat.mock.long.twophase.sample, assay == assays[aa]), aes_string(x = "Bserostatus", y = times[tt])) +
      geom_boxplot(aes(colour = Trt:Bserostatus), width = 0.6, lwd = 1.5) + 
      stat_boxplot(aes(colour = Trt:Bserostatus), geom = "errorbar", width = 0.45, lwd = 1.5) +
      guides(alpha = "none", fill = "none") +
      geom_jitter(data = subset(dat.mock.sample1, assay == assays[aa]), mapping = aes(colour = Trt:Bserostatus), width = 0.1, size = 5) +
      scale_y_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
      facet_wrap(vars(Trt), ncol = 1) + 
      theme_pubr(legend = "none") + 
      ylab(labels.axis[tt, aa]) + xlab("") + 
      scale_fill_manual(values = hvtn_col) +
      scale_color_manual(values = hvtn_col, labels = c("Placebo, Baseline Neg",
                                                       "Placebo, Baseline Pos",
                                                       "Vaccine, Baseline Neg",
                                                       "Vaccine, Baseline Pos")) +
      ggtitle(labels.title2[tt, aa]) +
      theme(plot.title = element_text(hjust = 0.5, size = 24),
            panel.border = element_rect(fill = NA),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_line(),
            axis.title = element_text(size = 24),
            axis.text = element_text(size = 28),
            strip.text = element_text(size = 20, face = "bold"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, face = "bold"))
    if (tt <= 3) {
      boxplot_list[[aa]] <- boxplot_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
        geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 8) 
    }
  }
  png(file = paste0(save.results.to, "boxplots_", times[tt], "_x_bstatus_", study.name, ".png"), height = 1720, width = 960)
  print(ggarrange(plotlist = boxplot_list, ncol = 2, nrow = 2, 
                  common.legend = TRUE, legend = "bottom",
                  align = "h")) 
  dev.off()
}

#=====================================================================
## Make seperate plots for Placebo and Vaccine arms

for (trt in 1:2) {
  for (tt in 1:5){
    boxplot_list <- vector("list", 4)
    for (aa in 1:4) {
      boxplot_list[[aa]] <- ggplot(subset(dat.mock.long.twophase.sample, assay == assays[aa] & as.numeric(Trt) == trt), aes_string(x = "Bserostatus", y = times[tt])) +
        geom_boxplot(aes(colour = Bserostatus), width = 0.6, lwd = 1.5) + 
        stat_boxplot(aes(colour = Bserostatus), geom = "errorbar", width = 0.45, lwd = 1.5) +
        guides(alpha = "none", fill = "none") +
        geom_jitter(data = subset(dat.mock.sample1, assay == assays[aa] & as.numeric(Trt) == trt), mapping = aes(colour = Bserostatus), width = 0.1, size = 5) +
        scale_y_continuous(limits = c(-2, 10), labels = label_math(10 ^ .x), breaks = seq(-2, 10, 2)) +
        theme_pubr(legend = "none") + 
        ylab(labels.axis[tt, aa]) + xlab("") + 
        scale_fill_manual(values = hvtn_col) +
        scale_color_manual(values = hvtn_col, labels = c("Baseline Neg", "Baseline Pos")) +
        ggtitle(labels.title2[tt, aa]) +
        theme(plot.title = element_text(hjust = 0.5, size = 24),
              panel.border = element_rect(fill = NA),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 28),
              strip.text = element_text(size = 20, face = "bold"),
              legend.title = element_blank(),
              legend.text = element_text(size = 24, face = "bold"))
      
      if (tt <= 3) {
        boxplot_list[[aa]] <- boxplot_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
          geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 8) 
      }
    }
    png(file = paste0(save.results.to, "boxplots_", times[tt], "_x_bstatus_", c("placebo_arm_", "vaccine_arm_")[trt],
                      study.name, ".png"), height = 1260, width = 1260)
    print(ggarrange(plotlist = boxplot_list, ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom",
                    align = "h")) 
    dev.off()
  }
}

#==========================================================================================================================================================================================================================
# The correlation of each antibody marker readout between Day 1 and Day 57, and between Day 1 and delta, 
# is examined within each of the baseline strata subgroups, and within each treatment arm and baseline serostratum.  
# Pairs plots/scatterplots will be used, annotated with baseline strata-adjusted Spearman rank correlations, implemented in the PResiduals R package \nocite{LiShepherd2012} available on CRAN. 
# For the correlation within each treatment arm and baseline serostratum, because PResiduals does not currently handle sampling weights, we will create a representative dataset through case deletion and apply PResiduals.
#==========================================================================================================================================================================================================================


#===================================
# plots within each baseline stratum
#===================================


# for (tt in 2:3) {
#   
#   for (bserostatus in 1:2) {
#     for (trt in 1:2) {    
#       for (aa in 1:4) {
#         
#         
#         
#         subdat <- subset(dat.mock.long, as.numeric(Bserostatus) == bserostatus & 
#                            as.numeric(Trt) == trt & as.numeric(assay) == aa)
#         lim <- c(-1, 6)
#         strata_corr <- sapply(1:max.stratum2, function(ss) {  ## baseline stratum adjusted Spearman rank correlation
#           corr <- try(round(cor(subdat[subdat$Bstratum2 == ss, "B"], subdat[subdat$Bstratum2 == ss, times2[tt]], 
#                                 use = "complete.obs", method = "spearman"), 2), silent = TRUE)
#           if (class(corr) == "try-error") corr <- NA
#           corr
#         })
#         png(filename = paste0(save.results.to, "stratum/scatterplots_", times[tt],"vB_", bstatus.labels.2[bserostatus], "_", trt.labels[trt], "_", assays[aa], ".png", sep = ""),
#             width = 720, height = 720)
#         pp <- ggscatter(x = "B", y = times2[tt], xlab = labels.axis[1, aa], ylab = labels.axis[tt, aa], xlim = lim, ylim = lim,
#                         data = subdat, title = assay.labels[aa]) %>%
#           facet("Bstratum2", panel.labs = list(Bstratum2 = paste0("Stratum ", 1:max.stratum,", cor: ", strata_corr))) +
#           geom_abline(slope = 1) +
#           theme(plot.title = element_text(hjust = 0.5, size = 20),
#                 panel.grid.minor.y = element_line(),
#                 panel.grid.major.y = element_line(),
#                 axis.title = element_text(size = 20),
#                 axis.text = element_text(size = 24),
#                 strip.text = element_text(size = 20, face = "bold"),
#                 legend.title = element_blank())
#         
#         pp <- pp + scale_x_continuous(breaks = c(0, 2, 4, 6), labels = label_math(10^.x)) + 
#           scale_y_continuous(breaks = c(0, 2, 4, 6), labels = label_math(10^.x))
#         
#         print(pp)
#         dev.off()
#       }
#     }
#   }
# }
# 
# 




#==============================================================================================
# Point estimates of Day 57 marker positive response rates for the vaccine arm by the baseline
# demographic subgroups and the baseline serostrata are provided, as well as pooled over base-
# line demographic strata. 95% confidence intervals about positive response rates are computed
# using the Wilson score method Agresti and Coull (1998). The point and 95% CI estimates
# include all of the data and use IPS weights.
#==============================================================================================
res_rate_under_sample <- dat.mock %>% 
  split(., list(.$Trt, .$Bserostatus)) %>%
  lapply(function(subdat) {
    propD29 <- sum(subdat$EventIndPrimaryD29) / nrow(subdat)
    ciD29 <- scoreci(sum(subdat$EventIndPrimaryD29), nrow(subdat), 0.95)
    propD57 <- sum(subdat$EventIndPrimaryD57) / nrow(subdat)
    ciD57 <- scoreci(sum(subdat$EventIndPrimaryD57), nrow(subdat), 0.95)
    data.frame(Trt = subdat$Trt[1], 
               Bserostatus = subdat$Bserostatus[1],
               propD29 = propD29, propD57 = propD57, 
               ci.lb.D29 = ciD29$conf.int[1], ci.ub.D29 = ciD29$conf.int[2],
               ci.lb.D57 = ciD57$conf.int[1], ci.ub.D57 = ciD57$conf.int[2])
  }) %>% bind_rows() %>%
  mutate(Trt = factor(Trt, levels = 0:1, labels = c("Placebo", "Vaccine")),
         Bserostatus = factor(Bserostatus, levels = 0:1,
                              labels = c("Baseline Neg", "Baseline Pos")))

png(filename = paste0(save.results.to, "barplot_responseRate_D29_trt_both_bserostatus_both_", study.name,".png"),
    width = 820, height = 820)
print(ggplot(data = res_rate_under_sample,
             aes(x = Trt:Bserostatus, y = propD29, color = Trt:Bserostatus, fill = Trt:Bserostatus, alpha = 0.1)) +
        geom_col(lwd = 1.4) +
        geom_errorbar(aes(ymin = ci.lb.D29, ymax = ci.ub.D29, color = Trt:Bserostatus), lwd = 1.4, width = 0.4) +
        xlab("") + ylab("D29 Positive Response Rate") + ylim(0, 0.1) +
        scale_x_discrete(labels = c("1A", "1B", "2A", "2B")) +
        scale_color_manual(values = hvtn_col) +
        scale_fill_manual(values = hvtn_col, labels = c("1A: Placebo, Baseline Neg",
                                                        "1B: Placebo, Baseline Pos",
                                                        "2A: Vaccine, Baseline Neg",
                                                        "2B: Vaccine, Baseline Pos")) +
        guides(alpha = "none", color = "none",
               fill = guide_legend(ncol = 2, byrow = TRUE)) + 
        theme_bw() + ggtitle("D29 Positive Response Rate by Treatment and Baseline Status") +
        theme(plot.title = element_text(hjust = 0.5, size = 24),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 20),
              legend.title = element_blank(),
              legend.text = element_text(size = 18),
              legend.position = "bottom"))
dev.off()

png(filename = paste0(save.results.to, "barplot_responseRate_D57_trt_both_bserostatus_both_", study.name,".png"),
    width = 820, height = 820)
print(ggplot(data = res_rate_under_sample,
             aes(x = Trt:Bserostatus, y = propD57, color = Trt:Bserostatus, fill = Trt:Bserostatus, alpha = 0.1)) +
        geom_col(lwd = 1.4) +
        geom_errorbar(aes(ymin = ci.lb.D57, ymax = ci.ub.D57, color = Trt:Bserostatus), lwd = 1.4, width = 0.4) +
        xlab("") + ylab("D29 Positive Response Rate") + ylim(0, 0.1) +
        scale_x_discrete(labels = c("1A", "1B", "2A", "2B")) +
        scale_color_manual(values = hvtn_col) +
        scale_fill_manual(values = hvtn_col, labels = c("1A: Placebo, Baseline Neg",
                                                        "1B: Placebo, Baseline Pos",
                                                        "2A: Vaccine, Baseline Neg",
                                                        "2B: Vaccine, Baseline Pos")) +
        guides(alpha = "none", color = "none",
               fill = guide_legend(ncol = 2, byrow = TRUE)) + 
        theme_bw() + ggtitle("D57 Positive Response Rate by Treatment and Baseline Status") +
        theme(plot.title = element_text(hjust = 0.5, size = 24),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 20),
              legend.title = element_blank(),
              legend.text = element_text(size = 18),
              legend.position = "bottom"))
dev.off()
#==========================================================================================================================================================================================================================
# The correlation of each antibody marker readout between Day 1 and Day 57, and between Day 1 and delta, 
# is examined within each of the baseline strata subgroups, and within each treatment arm and baseline serostratum.  
# Pairs plots/scatterplots will be used, annotated with baseline strata-adjusted Spearman rank correlations, implemented in the PResiduals R package \nocite{LiShepherd2012} available on CRAN. 
# For the correlation within each treatment arm and baseline serostratum, because PResiduals does not currently handle sampling weights, we will create a representative dataset through case deletion and apply PResiduals.
#==========================================================================================================================================================================================================================


#===================================
# plots within each baseline stratum
#===================================


# for (tt in 2:5) {
#   
#   for (bserostatus in 1:2) {
#     for (trt in 1:2) {    
#       for (aa in 1:4) {
#         
#         
#         
#         subdat <- subset(dat.mock.long, as.numeric(Bserostatus) == bserostatus & 
#                            as.numeric(Trt) == trt & as.numeric(assay) == aa)
#         lim <- c(-1, 6)
#         strata_corr <- sapply(1:8, function(ss) {  ## baseline stratum adjusted Spearman rank correlation
#           round(cor(subdat[subdat$Bstratum == ss, "B"], subdat[subdat$Bstratum == ss, times[tt]], 
#                     use = "complete.obs", method = "spearman"), 2)
#         })
#         png(filename = paste0(save.results.to, "stratum/scatterplots_", times[tt],"vB_", bstatus.labels.2[bserostatus], "_", trt.labels[trt], "_", assays[aa], ".png", sep = ""),
#             width = 960, height = 960)
#         pp <- ggscatter(x = "B", y = times[tt], xlab = labels.axis[1, aa], ylab = labels.axis[tt, aa], xlim = lim, ylim = lim,
#                         data = subdat,  title = assay.labels[aa]) %>%
#           facet("Bstratum", panel.labs = list(Bstratum = paste0("Stratum ", 1:8, ", cor: ", strata_corr))) +
#           geom_abline(slope = 1) +
#           theme(plot.title = element_text(hjust = 0.5, size = 24),
#                 panel.border = element_rect(fill = NA),
#                 panel.grid.minor.y = element_line(),
#                 panel.grid.major.y = element_line(),
#                 axis.title = element_text(size = 20),
#                 axis.text = element_text(size = 28),
#                 strip.text = element_text(size = 20, face = "bold"),
#                 legend.title = element_blank())
#         
#         pp <- pp + scale_x_continuous(breaks = c(0, 2, 4, 6), labels = label_math(10^.x)) + 
#           scale_y_continuous(breaks = c(0, 2, 4, 6), labels = label_math(10^.x))
#         
#         print(pp)
#         dev.off()
#       }
#     }
#   }
# }







##===========================================================================================================================================================
## Boxplots restricting to vaccine recipients who are non-Cases and focusing on the one of the Day 57 marker readouts not subtracting for baseline.  
## That is, 4 of these figures, one each for each the 4 Day 57 antibody markers not subtracting for baseline.
## The subgroups are: Age < 65, Age >= 65, High Risk and Not High Risk (for the pre-existing high risk conditions variable), Minority and non-Minority.
## Objective of figure: describe if the antibody markers in vaccine recipients tends to differ over major subgroups of interest.  The plot restricts to 
## non-Cases as a simple way to approximately show the distribution for the whole cohort.
##===========================================================================================================================================================


# 
# dat.mock.sample2 <- dat.mock.long[, c("Trt", "Bserostatus", "B", "Day57", "Delta", "assay", "Bstratum2",
#                                       "EventInd")] %>%
#   filter(., complete.cases(.)) %>% 
#   split(., list(.$Trt, .$Bserostatus, .$Bstratum2, .$EventInd)) %>%
#   lapply(., function(x) {
#     if(nrow(x) <= 30) {
#       return(x)
#     } else {
#       return(x[sample(1:nrow(x), size = 30),])
#     }}) %>% bind_rows
# 
# 
# for (bstatus in 1:2) {
#   for (tt in 2:3){
#     boxplot_list <- vector("list", 4)
#     for (aa in 1:4) {
#       boxplot_list[[aa]] <- ggplot(subset(dat.mock.long, assay == assays[aa] & Trt == "Vaccine" & EventInd == 0 &
#                                             Bserostatus == bstatus.labels[bstatus]), aes_string(x = "Bstratum2", y = times2[tt])) +
#         geom_boxplot(aes(colour = Bstratum2), width = 0.6, lwd = 1.5) + 
#         stat_boxplot(aes(colour = Bstratum2), geom = "errorbar", width = 0.45, lwd = 1.5) +
#         geom_jitter(data = filter(dat.mock.sample2, assay == assays[aa] & Trt == "Vaccine" & EventInd == 0 & 
#                                     Bserostatus == bstatus.labels[bstatus]), mapping = aes(colour = Bstratum2), width = 0.1, size = 3) +
#         theme_pubr(legend = "none") + 
#         ylab(labels.axis[tt, aa]) + xlab("") + ggtitle(labels.title[tt, aa]) +
#         guides(alpha = "none") +
#         scale_color_manual(values = hvtn_col,
#                            labels = bstrata.labels) +
#         scale_y_continuous(limits = c(-1, 6), labels = label_math(10^.x)) +
#         theme(plot.title = element_text(hjust = 0.5, size = 20),
#               panel.border = element_rect(fill = NA),
#               panel.grid.minor.y = element_line(),
#               panel.grid.major.y = element_line(),
#               axis.title = element_text(size = 20),
#               axis.text.x = element_blank(),
#               axis.text.y = element_text(size = 28),
#               strip.text = element_text(size = 20, face = "bold"),
#               legend.title = element_blank(),
#               legend.text = element_text(size = 20))
#     }
#     png(file = paste0(save.results.to, "boxplots_", times2[tt], "_trt_vaccine_controls_stratified_noncase_",
#                       bstatus.labels.2[bstatus], ".png"), height = 960, width = 1560)
#     print(ggarrange(plotlist = boxplot_list, ncol = 2, nrow = 2, 
#                     common.legend = TRUE, legend = "right", 
#                     align = "h")  + 
#             theme(plot.title = element_text(hjust = 0.4, size = 24),
#                   legend.title = element_blank(),
#                   legend.text = element_text(size = 14)))
#     dev.off()
#   }
# }

##===========================================================================
## Box plots comparing the day 57 assay readouts between case/non-case and 
## vaccine-placebo for four different assays
##===========================================================================
## make another subsample datasets such that the jitter plot for each subgroup in each panel contains no more
## than 50 data points
dat.mock.sample3 <-  dat.mock.long.twophase.sample[, c("Trt", "Bserostatus", "B", "Day29", "Day57", "Delta29overB", "Delta57overB", "assay", "EventD57")] %>%
  filter(., .$Trt == "Vaccine") %>% 
  split(., list(.$Bserostatus, .$assay, .$EventD57)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows


dat.mock.sample4 <-  dat.mock.long[, c("Trt", "Bserostatus", "B", "Day29", "Day57", "Delta29overB", "Delta57overB", "assay", "EventD57")] %>%
  filter(., .$Trt == "Vaccine") %>% 
  split(., list(.$Bserostatus, .$assay, .$EventD57)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows
##===============================================================================================================
## Box plots, overlayed with boxplots and jittered sample points, for the distribution of day 57 assay readouts
## versus the case-control status among the vaccine recipients, stratified by the baseline sero-status
##===============================================================================================================


for (bstatus in 1:2) {
  for (tt in 2:5){
    boxplot_list <- vector("list", 4)
    for (aa in 1:4) {
      boxplot_list[[aa]] <-  ggplot(subset(dat.mock.long, assay == assays[aa] & Trt == "Vaccine" & 
                                             Bserostatus == bstatus.labels[bstatus]), aes_string(x = "EventD57", y = times[tt])) +
        geom_boxplot(aes(colour = EventD57), width = 0.6, lwd = 1.5) + 
        stat_boxplot(geom = "errorbar", aes(colour = EventD57), width = 0.45, lwd = 1.5) +
        geom_jitter(data = filter(dat.mock.sample4, assay == assays[aa] & Trt == "Vaccine" & 
                                    Bserostatus == bstatus.labels[bstatus]), mapping = aes(colour = EventD57), width = 0.1, 
                    size = 5, show.legend = FALSE) +
        theme_pubr(legend = "none") + 
        guides(alpha = "none") +
        ylab(labels.axis[tt, aa]) + xlab("") + ggtitle(labels.title2[tt, aa]) +
        scale_color_manual(values = hvtn_col,
                           labels = c("D57 Non-Case", "D57 Case")) +
        scale_y_continuous(limits = c(-2, 6), labels = label_math(10^.x), breaks = seq(-2, 10, 2)) +
        theme(plot.title = element_text(hjust = 0.5, size = 24),
              panel.border = element_rect(fill = NA),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 20),
              axis.text.x = element_text(size = 28),
              axis.text.y = element_text(size = 28),
              strip.text = element_text(size = 20, face = "bold"),
              legend.title = element_blank(),
              legend.text = element_text(size = 24),
              legend.spacing = unit(1, "cm"))
      if (tt <= 3) {
        boxplot_list[[aa]] <- boxplot_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
          geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 8) 
      }
    }
    png(file = paste0(save.results.to, "boxplots_", times[tt], "_trt_vaccine_x_cc_", bstatus.labels.2[bstatus], 
                      "_", study.name, ".png"), 
        height = 1260, width = 1260)
    print(ggarrange(plotlist = boxplot_list, 
                    ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom",
                    align = "h") + 
            theme(plot.title = element_text(hjust = 0.5, size = 24)))
    
    dev.off()
  }
}









##===============================================================================================================
## Box plots, overplayed with boxplots and jittered sample points, for the distribution of baseline, day 57, or 
## fold-rise in assay readouts versus the case-control status among the vaccine recipients, 
## stratified by the baseline serostatus, age group and risk group
##===============================================================================================================

dat.mock.sample5 <-  dat.mock.long[, c("Trt", "Bserostatus", "B", "Day29", "Day57", "Delta29overB", 
                                       "Delta57overB", "assay", "EventD57", "demo_lab")] %>%
  filter(., .$Trt == "Vaccine") %>% 
  split(., list(.$Bserostatus, .$assay, .$EventD57, .$demo_lab)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows


for (bstatus in 1:2) {
  for (tt in 2:5){
    for (aa in 1:4) {
      boxplots <-  ggplot(subset(dat.mock.long, assay == assays[aa] & Trt == "Vaccine" & 
                                   Bserostatus == bstatus.labels[bstatus]), aes_string(x = "EventD57", y = times[tt])) +
        geom_boxplot(aes(colour = EventD57), width = 0.6, lwd = 1.5) + 
        stat_boxplot(geom = "errorbar", aes(colour = EventD57), width = 0.45, lwd = 1.5) +
        geom_jitter(data = filter(dat.mock.sample5, assay == assays[aa] & 
                                    Bserostatus == bstatus.labels[bstatus]), mapping = aes(colour = EventD57), width = 0.1, 
                    size = 5, show.legend = FALSE) +
        theme_pubr() + 
        guides(alpha = "none") +
        facet_wrap(~ demo_lab) +
        ylab(labels.axis[tt, aa]) + xlab("") + ggtitle(labels.title2[tt, aa]) +
        scale_color_manual(values = hvtn_col,
                           labels = c("D57 Non-Case", "D57 Case")) +
        scale_y_continuous(limits = c(-2, 10), labels = label_math(10^.x), breaks = seq(-2, 10, 2)) +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
              panel.border = element_rect(fill = NA),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 20),
              axis.text.x = element_text(size = 28),
              axis.text.y = element_text(size = 28),
              strip.text = element_text(size = 20, face = "bold"),
              legend.position = "bottom",
              legend.title = element_blank(),
              legend.text = element_text(size = 24),
              legend.spacing = unit(1, "cm"))
      
      if (tt <= 3) {
        boxplots <- boxplots + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
          geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 8) 
      }
      
      png(file = paste0(save.results.to, "boxplots_", times[tt],"_trt_vaccine_x_cc_", 
                        bstatus.labels.2[bstatus],"_", assays[aa], "_", study.name, ".png"), 
          height = 960, width = 960)
      print(boxplots)
      dev.off()
    }
    
  }
}

