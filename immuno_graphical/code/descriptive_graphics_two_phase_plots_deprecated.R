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


