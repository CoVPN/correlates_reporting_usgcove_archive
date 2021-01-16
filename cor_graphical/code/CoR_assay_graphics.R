source("descriptive_graphics_data_preprocess.R")
source("ggally_cor_resample.R") 


#=========================================================================================================================
# Reverse empirical cdf (rcdf) plots for the Baseline/Day 57/Baseline-subtracted Day 57assay readouts, 
# stratified by treatment group and event status, in baseline negative or positive subjects
# We made four ggplot objects, each for one assay, and combine them with ggarrange
#=========================================================================================================================
for (eventDay in c(29, 57)) {
  for (bstatus in 1:2) {
    for (tt in 1:5){
      rcdf_list <- vector("list", 4)
      for (aa in 1:4) {
        rcdf_list[[aa]] <- ggplot(subset(dat.long.twophase.sample, assay == assays[aa] & Bserostatus == bstatus.labels[bstatus]), 
                                  aes_string(x = times[tt], colour = paste0("Trt:EventLabelD", eventDay), 
                                             linetype = paste0("EventLabelD", eventDay), weight = "wt")) +
          geom_line(aes(y = 1 - ..y..), stat='ecdf', lwd = 1) +  theme_pubr(legend = "none") +
          ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) +
          scale_x_continuous(labels = label_math(10^.x), limits = c(-2, 10), breaks = seq(-2, 10, 2)) +
          scale_color_manual(values = hvtn_col) +
          guides(linetype = "none",
                 color = guide_legend(nrow = 2, byrow = TRUE)) +
          ggtitle(labels.title2[tt, aa]) +
          theme(plot.title = element_text(hjust = 0.5, size = 10),
                legend.title = element_blank(),
                legend.text = element_text(size = 10, face = "bold"),
                panel.grid.minor.y = element_line(),
                panel.grid.major.y = element_line(),
                axis.title = element_text(size = 9),
                axis.text = element_text(size = 10))
      }
      
      ggsave(ggarrange(plotlist = rcdf_list, ncol = 2, nrow = 2, 
                       common.legend = TRUE, legend = "bottom",
                       align = "h"),
             filename = paste0(save.results.to, "Marker_RCDF_", times[tt], 
                               "_trt_by_caseD", eventDay,"_bstatus_", bstatus.labels.2[bstatus],
                               "_", study.name,".png"),
             height = 7, width = 6.5)
      
    }
  }
}


##===========================================================================
## Box plots comparing the day 57 assay readouts between case/non-case and 
## vaccine-placebo for four different assays
##===========================================================================
## make another subsample datasets such that the jitter plot for each subgroup in each panel contains no more
## than 50 data points
dat.sample3 <-  dat.long.twophase.sample[, c("Trt", "Bserostatus", times, "assay", "EventD57")] %>%
  filter(., .$Trt == "Vaccine") %>% 
  split(., list(.$Bserostatus, .$assay, .$EventD57)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows


dat.sample4 <-  dat.long[, c("Trt", "Bserostatus", times, "assay", "EventD57")] %>%
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
      boxplot_list[[aa]] <-  ggplot(subset(dat.long, assay == assays[aa] & Trt == "Vaccine" & 
                                             Bserostatus == bstatus.labels[bstatus]), aes_string(x = "EventD57", y = times[tt])) +
        geom_boxplot(aes(colour = EventD57), width = 0.6, lwd = 1) + 
        stat_boxplot(geom = "errorbar", aes(colour = EventD57), width = 0.45, lwd = 1) +
        geom_jitter(data = filter(dat.sample4, assay == assays[aa] & Trt == "Vaccine" & 
                                    Bserostatus == bstatus.labels[bstatus]), mapping = aes(colour = EventD57), width = 0.1, 
                    size = 1.4, show.legend = FALSE) +
        theme_pubr(legend = "none") + 
        guides(alpha = "none") +
        ylab(labels.axis[tt, aa]) + xlab("") + ggtitle(labels.title2[tt, aa]) +
        scale_color_manual(values = hvtn_col,
                           labels = c("D57 Non-Case", "D57 Case")) +
        scale_y_continuous(limits = c(-2, 10), labels = label_math(10^.x), breaks = seq(-2, 10, 2)) +
        theme(plot.title = element_text(hjust = 0.5, size = 10),
              panel.border = element_rect(fill = NA),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 9),
              axis.text.x = element_text(size = 10),
              axis.text.y = element_text(size = 10),
              strip.text = element_text(size = 10, face = "bold"),
              legend.title = element_blank())
      if (tt <= 3) {
        boxplot_list[[aa]] <- boxplot_list[[aa]] + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
          geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 5) 
      }
    }

    ggsave(ggarrange(plotlist = boxplot_list, 
                    ncol = 2, nrow = 2, 
                    common.legend = TRUE, legend = "bottom",
                    align = "h") + 
            theme(plot.title = element_text(hjust = 0.5, size = 10)),
           filename = paste0(save.results.to, "boxplots_", times[tt], "_trt_vaccine_x_cc_", bstatus.labels.2[bstatus], 
                             "_", study.name, ".png"),
           height = 7, width = 6.5)

  }
}









##===============================================================================================================
## Box plots, overplayed with boxplots and jittered sample points, for the distribution of baseline, day 57, or 
## fold-rise in assay readouts versus the case-control status among the vaccine recipients, 
## stratified by the baseline serostatus, age group and risk group
##===============================================================================================================

dat.sample5 <-  dat.long[, c("Trt", "Bserostatus", times, "assay", "EventD57", "demo_lab")] %>%
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
      boxplots <-  ggplot(subset(dat.long, assay == assays[aa] & Trt == "Vaccine" & 
                                   Bserostatus == bstatus.labels[bstatus]), aes_string(x = "EventD57", y = times[tt])) +
        geom_boxplot(aes(colour = EventD57), width = 0.6, lwd = 1) + 
        stat_boxplot(geom = "errorbar", aes(colour = EventD57), width = 0.45, lwd = 1) +
        geom_jitter(data = filter(dat.sample5, assay == assays[aa] & 
                                    Bserostatus == bstatus.labels[bstatus]), mapping = aes(colour = EventD57), width = 0.1, 
                    size = 1.4, show.legend = FALSE) +
        theme_pubr() + 
        guides(alpha = "none") +
        facet_wrap(~ demo_lab) +
        ylab(labels.axis[tt, aa]) + xlab("") + ggtitle(labels.title2[tt, aa]) +
        scale_color_manual(values = hvtn_col,
                           labels = c("D57 Non-Case", "D57 Case")) +
        scale_y_continuous(limits = c(-2, 10), labels = label_math(10^.x), breaks = seq(-2, 10, 2)) +
        theme(plot.title = element_text(hjust = 0.5, size = 10),
              panel.border = element_rect(fill = NA),
              panel.grid.minor.y = element_line(),
              panel.grid.major.y = element_line(),
              axis.title = element_text(size = 9),
              axis.text.x = element_text(size = 10),
              axis.text.y = element_text(size = 10),
              strip.text = element_text(size = 10, face = "bold"),
              legend.position = "bottom",
              legend.title = element_blank(),
              legend.text = element_text(size = 10))
      
      if (tt <= 3) {
        boxplots <- boxplots + geom_hline(yintercept = LLOQ[aa], linetype = 2, color = "black", lwd = 1) +
          geom_text(x = 0.7, vjust = "right", y = LLOQ[aa] - 0.5 , label = "LLOQ", size = 5) 
      }
      

      ggsave(boxplots,
             filename = paste0(save.results.to, "boxplots_", times[tt],"_trt_vaccine_x_cc_", 
                               bstatus.labels.2[bstatus],"_", assays[aa], "_", study.name, ".png"),
             height = 6.5,
             width = 6.5)
    }
  }
}

