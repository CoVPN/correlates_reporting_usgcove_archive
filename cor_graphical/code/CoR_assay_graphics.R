source("descriptive_graphics_data_preprocess.R")
source("ggally_cor_resample.R") 


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

