#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------
source(here::here("code", "params.R"))
source(here::here("code", "cor_data_preprocess.R"))
source(here::here("code", "covid_corr_plot_functions.R"))
library(ggpubr)
library(scales)
library(ggplot2)
#=========================================================================================================================
# Reverse empirical cdf (rcdf) plots for the Baseline/Day 57/Baseline-subtracted Day 57assay readouts, 
# stratified by treatment group and event status, in baseline negative or positive subjects
# We made four ggplot objects, each for one assay, and combine them with ggarrange
#=========================================================================================================================
wts <- c("wt", "wt.2", "wt", "wt.2", "wt")
for (bstatus in 1:2) {
  for (tt in 2:5){
    subdat <- subset(dat.long.cor.subset, Bserostatus == bstatus.labels[bstatus])
    subdat$TrtEvent <- paste(as.character(subdat$Trt), as.character(subdat$cohort_event), ", ")
    if (tt %in% c(3, 5)) {  ## day 57 analysis don't include intercurrent cases
      subdat <- subdat %>% filter(cohort_event != "Intercurrent Cases")
    }
    rcdf_list <- vector("list", 4)
    for (aa in 1:4) {
      rcdf_list[[aa]] <- ggplot(subset(subdat, assay == assays[aa]), 
                                aes_string(x = times[tt], colour = "TrtEvent", 
                                           linetype = "TrtEvent", weight = wts[tt])) +
        geom_step(aes(y = 1 - ..y..), stat = "ecdf", lwd = 1) +
        theme_pubr(legend = "none") +
        ylab("Reverse ECDF") + xlab(labels.axis[tt, aa]) +
        scale_x_continuous(labels = label_math(10^.x), limits = c(-2, 10), breaks = seq(-2, 10, 2)) +
        scale_color_manual(values = c("#1749FF", "#D92321", "#0AB7C9", "#FF6F1B", "#810094", "#378252", "#FF5EBF", "#3700A5", "#8F8F8F", "#787873")) +
        guides(linetype = "none",
               color = guide_legend(nrow = 3, byrow = TRUE)) +
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
           filename = paste0(save.results.to, "/Marker_RCDF_", times[tt], 
                             "_trt_by_event_status_bstatus_", bstatus.labels.2[bstatus], "_", 
                             study.name,".png"),
           height = 7, width = 6.5)
    
    
  }
}


##===========================================================================
## Box plots comparing the day 57 assay readouts between case/non-case and 
## vaccine-placebo for four different assays
##===========================================================================
## make another subsample datasets such that the jitter plot for each subgroup in each panel contains no more
## than 50 data points
set.seed(12345)
dat.sample3 <-  dat.long.cor.subset[, c("Trt", times, "assay", "cohort_event", "Bserostatus")] %>%
  filter(., .$Trt == "Vaccine") %>% 
  split(., list(.$assay, .$cohort_event, .$Bserostatus)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows

set.seed(12345)
dat.sample4 <-  dat.long.cor.subset[, c("Trt", times, "assay", "cohort_event", "Bserostatus")] %>%
  filter(., .$Trt == "Vaccine") %>% 
  split(., list(.$assay, .$cohort_event, .$Bserostatus)) %>%
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
    subdat <- subset(dat.long.cor.subset, Bserostatus == bstatus.labels[bstatus])
    subdat_jitter <- subset(dat.sample4, Bserostatus == bstatus.labels[bstatus])
    if (tt %in% c(3, 5)) {
      subdat <- subdat %>% filter(cohort_event != "Intercurrent Cases")
      subdat_jitter <- subdat_jitter %>% filter(cohort_event != "Intercurrent Cases")
    }
    boxplot_list <- vector("list", 4)
    for (aa in 1:4) {
      boxplot_list[[aa]] <-  ggplot(subset(subdat, assay == assays[aa] & Trt == "Vaccine"), 
                                    aes_string(x = "cohort_event", y = times[tt])) +
        geom_boxplot(aes(colour = cohort_event), width = 0.6, lwd = 1) + 
        stat_boxplot(geom = "errorbar", aes(colour = cohort_event), width = 0.45, lwd = 1) +
        geom_jitter(data = filter(subdat_jitter, assay == assays[aa] & Trt == "Vaccine"), mapping = aes(colour = cohort_event), width = 0.1, 
                    size = 1.4, show.legend = FALSE) +
        theme_pubr(legend = "none") + 
        guides(alpha = "none") +
        ylab(labels.axis[tt, aa]) + xlab("") + ggtitle(labels.title2[tt, aa]) +
        scale_color_manual(values = c("#1749FF", "#D92321", "#0AB7C9", "#FF6F1B", "#810094", "#378252", "#FF5EBF", "#3700A5", "#8F8F8F", "#787873")) +
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
        boxplot_list[[aa]] <- boxplot_list[[aa]] + geom_hline(yintercept = log10(llods[aa]), linetype = 2, color = "black", lwd = 1) +
          geom_text(x = 0.7, vjust = "right", y = log10(llods[aa]) - 0.5 , label = "LLOD", size = 5) 
      }
    }
    
    ggsave(ggarrange(plotlist = boxplot_list, 
                     ncol = 2, nrow = 2, 
                     common.legend = TRUE, legend = "bottom",
                     align = "h") + 
             theme(plot.title = element_text(hjust = 0.5, size = 10)),
           filename = paste0(save.results.to, "/boxplots_", times[tt], "_trt_vaccine_x_cc_",
                             bstatus.labels.2[bstatus], "_",
                             study.name, ".png"),
           height = 9, width = 8)
    
  }
}









##===============================================================================================================
## Box plots, overlayed with boxplots and jittered sample points, for the distribution of baseline, day 57, or 
## fold-rise in assay readouts versus the case-control status among the vaccine recipients, 
## stratified by the baseline serostatus, age group and risk group
##===============================================================================================================
set.seed(12345)
dat.sample5 <-  dat.long.cor.subset[, c("Trt", times, "assay", "cohort_event", "demo_lab", "Bserostatus")] %>%
  filter(., .$Trt == "Vaccine") %>% 
  split(., list(.$assay, .$cohort_event, .$demo_lab, .$Bserostatus)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 30) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 30),])
    }}) %>% bind_rows


for (bstatus in 1:2) {
  for (tt in 2:5){
    subdat <- subset(dat.long.cor.subset, Bserostatus == bstatus.labels[bstatus])
    subdat_jitter <- subset(dat.sample4, Bserostatus == bstatus.labels[bstatus])
    if (tt %in% c(3, 5)) {
      subdat <- subdat %>% filter(cohort_event != "Intercurrent Cases")
      subdat_jitter <- subdat_jitter %>% filter(cohort_event != "Intercurrent Cases")
    }
    for (aa in 1:4) {
      boxplots <-  ggplot(subset(subdat, assay == assays[aa] & Trt == "Vaccine"), aes_string(x = "cohort_event", y = times[tt])) +
        geom_boxplot(aes(colour = cohort_event), width = 0.6, lwd = 1) + 
        stat_boxplot(geom = "errorbar", aes(colour = cohort_event), width = 0.45, lwd = 1) +
        geom_jitter(data = filter(subdat_jitter, assay == assays[aa] & Trt == "Vaccine"), mapping = aes(colour = cohort_event), width = 0.1, 
                    size = 1.4, show.legend = FALSE) +
        theme_pubr() + 
        guides(alpha = "none", color = guide_legend(nrow = 2, byrow = TRUE)) +
        facet_wrap(~ demo_lab) +
        ylab(labels.axis[tt, aa]) + xlab("") + ggtitle(labels.title2[tt, aa]) +
        scale_color_manual(values = c("#1749FF", "#D92321", "#0AB7C9", "#FF6F1B", "#810094", "#378252", "#FF5EBF", "#3700A5", "#8F8F8F", "#787873")) +
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
        boxplots <- boxplots + geom_hline(yintercept = log10(llods[aa]), linetype = 2, color = "black", lwd = 1) +
          geom_text(x = 0.7, vjust = "right", y = log10(llods[aa]) - 0.5 , label = "LLOD", size = 5) 
      }
      
      
      ggsave(boxplots,
             filename = paste0(save.results.to, "/boxplots_", times[tt],"_trt_vaccine_x_cc_",
                               bstatus.labels.2[bstatus], "_", 
                               assays[aa], "_", study.name, ".png"),
             height = 9,
             width = 8)
    }
  }
}

#-----------------------------------------------
# Spaghetti PLOTS
#-----------------------------------------------
# - Spaghetti plots of antibody marker change over time
#-----------------------------------------------

## in each baseline serostatus group, randomly select 10 placebo recipients and 20 vaccine recipients
set.seed(12345)
var_names <- expand.grid(times = c("B", "Day29", "Day57"),
                         assays = assays) %>%
  mutate(var_names = paste0(times, assays)) %>%
  .[, "var_names"]

spaghetti_ptid <- dat.cor.subset[, c("Ptid", "Bserostatus", "Trt", var_names, "cohort_event")] %>%
  filter(., complete.cases(.), Trt == 1) %>%
  transmute(Bserostatus = Bserostatus,
            cohort_event = cohort_event,
            Ptid = Ptid) %>%
  split(., list(.$Bserostatus, .$cohort_event)) %>%
  lapply(function(xx) {
    if (nrow(xx) <= 15) {
      return(xx$Ptid)
    } else {
      return(sample(xx$Ptid, 15))
    }
  }) %>% unlist %>% as.numeric

spaghetti_dat <- dat.long.cor.subset[, c("Ptid", "Bserostatus", "cohort_event", 
                                         "B", "Day29", "Day57", "assay")] %>%
  filter(Ptid %in% spaghetti_ptid) %>%
  pivot_longer(cols = c("B", "Day29", "Day57"),
               names_to = "time") %>%
  mutate(assay_label = factor(assay, levels = assays, labels = labels.assays.short[assays]),
         time_label = factor(time, levels = c("B", "Day29", "Day57"),
                             labels = c("D1", "D29", "D57"))) %>%
  as.data.frame

for (bstatus in 1:2) {
  subdat <- subset(spaghetti_dat, Bserostatus == bstatus.labels[bstatus])
  covid_corr_spaghetti_facets(plot_dat = subdat,
                              x = "time_label",
                              y = "value",
                              id = "Ptid",
                              color = "cohort_event",
                              facet_by = "assay_label",
                              plot_title = paste0(
                                "Baseline ",
                                c("Negative", "Positive")[bstatus],
                                " PP Vaccine group"
                              ),
                              filename = paste0(
                                save.results.to, "/spaghetti_plot_trt_",
                                bstatus.labels.2[bstatus], "_",
                                study.name, ".png"
                              ),
                              height = 6, width = 5)
}
