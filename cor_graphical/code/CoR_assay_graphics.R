#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------
source(here::here("code", "params.R"))
source(here::here("code", "covid_corr_plot_functions.R"))
library(ggpubr)
library(scales)
library(ggplot2)
library(tidyr)

dat.long.cor.subset <- readRDS(here(
  "data_clean",
  "long_cor_data.rds"
))

dat.cor.subset <- readRDS(here(
  "data_clean",
  "cor_data.rds"
))

#=========================================================================================================================
# Reverse empirical cdf (rcdf) plots for the Baseline/Day 57/Baseline-subtracted Day 57assay readouts, 
# stratified by treatment group and event status, in baseline negative or positive subjects
# We made four ggplot objects, each for one assay, and combine them with ggarrange
#=========================================================================================================================
wts <- c("wt.D57", "wt.D29", "wt.D57", "wt.D29", "wt.D57")
tps <- c("Day29", "Day57", "Delta29overB", "Delta57overB")
for (tp in tps[tps %in% times]){
  dat.long.cor.subset$TrtEvent <- paste(as.character(dat.long.cor.subset$Trt), as.character(dat.long.cor.subset$cohort_event),sep = ", ")
  if (tp %in% c("Day57", "Delta57overB")) {  ## day 57 analysis don't include intercurrent cases
    dat.long.cor.subset <- dat.long.cor.subset %>% filter(cohort_event != "Intercurrent Cases" & ph2.D57==1)
  }
  rcdf_list <- vector("list", 4)
  for (aa in seq_along(assays)) {
    rcdf_list[[aa]] <- ggplot(subset(dat.long.cor.subset, assay == assays[aa]), 
                              aes_string(
                                x = tp, 
                                colour = "TrtEvent", 
                                linetype = "cohort_event",
                                weight = ifelse(
                                  tp %in% c("Day29", "Delta29overB"),
                                  "wt.D29", 
                                  "wt.D57")
                              )
    ) +
      geom_step(aes(y = 1 - ..y..), stat = "ecdf", lwd = 1) +
      theme_pubr(legend = "none") +
      ylab("Reverse ECDF") + xlab(labels.axis[tp, aa]) +
      scale_x_continuous(labels = label_math(10^.x), limits = c(-2, 6), breaks = seq(-2, 6, 2)) +
      scale_color_manual(values = c("#1749FF", "#D92321", "#0AB7C9", "#FF6F1B", "#810094", "#378252", "#FF5EBF", "#3700A5", "#8F8F8F", "#787873")) +
      guides(linetype = "none",
             color = guide_legend(nrow = 3, byrow = TRUE)) +
      ggtitle(labels.title2[tp, aa]) +
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            legend.title = element_blank(),
            legend.text = element_text(size = 14),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_line(),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 14))
  }
  
  ggsave(ggarrange(plotlist = rcdf_list, ncol = 2, nrow = 2, 
                   common.legend = TRUE, legend = "bottom",
                   align = "h"),
         filename = paste0(save.results.to, "/Marker_RCDF_", tp, 
                           "_trt_by_event_status_bstatus_", 
                           study_name,".png"),
         height = 7, width = 6.5)
  
  
}



##===========================================================================
## Box plots comparing the day 57 assay readouts between case/non-case and 
## vaccine-placebo for four different assays
##===========================================================================
## make another subsample datasets such that the jitter plot for each subgroup in each panel contains no more
## than 50 data points
set.seed(12345)
dat.sample3 <- dat.long.cor.subset %>%
  filter(., .$Trt == "Vaccine") %>% 
  split(., list(.$assay, .$cohort_event)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 100) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 100),])
    }}) %>% bind_rows

set.seed(12345)
dat.sample4 <-  dat.long.cor.subset %>%
  filter(., .$Trt == "Vaccine") %>% 
  split(., list(.$assay, .$cohort_event)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 100) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 100),])
    }}) %>% bind_rows
##===============================================================================================================
## Box plots, overlayed with boxplots and jittered sample points, for the distribution of day 57 assay readouts
## versus the case-control status among the vaccine recipients, stratified by the baseline sero-status
##===============================================================================================================



for (tp in tps[tps %in% times]){
  subdat <- dat.long.cor.subset
  subdat_jitter <- dat.sample4
  if (tp %in% c("Day57", "Delta57overB")) {
    subdat <- dat.long.cor.subset %>% filter(cohort_event != "Intercurrent Cases" & ph2.D57==1)
    subdat_jitter <- subdat_jitter %>% filter(cohort_event != "Intercurrent Cases" & ph2.D57==1)
  }
  boxplot_list <- vector("list", 4)
  for (aa in seq_along(assays)) {
    boxplot_list[[aa]] <-  ggplot(subset(subdat, assay == assays[aa] & Trt == "Vaccine"), 
                                  aes_string(x = "cohort_event", y = tp)) +
      geom_boxplot(aes(colour = cohort_event), width = 0.6, lwd = 1) + 
      stat_boxplot(geom = "errorbar", aes(colour = cohort_event), width = 0.45, lwd = 1) +
      geom_jitter(data = filter(subdat_jitter, assay == assays[aa] & Trt == "Vaccine"), 
                  mapping = aes(colour = cohort_event), width = 0.1, 
                  size = 1.4, show.legend = FALSE) +
      theme_pubr(legend = "none") + 
      guides(alpha = "none") +
      ylab(labels.axis[tp, aa]) + xlab("") + ggtitle(labels.title2[tp, aa]) +
      scale_color_manual(values = c("#1749FF", "#D92321", "#0AB7C9", "#FF6F1B", "#810094", "#378252", 
                                    "#FF5EBF", "#3700A5", "#8F8F8F", "#787873")) +
      scale_y_continuous(limits = c(-2, 6), labels = label_math(10^.x), breaks = seq(-2, 6, 2)) +
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            panel.border = element_rect(fill = NA),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_line(),
            axis.title = element_text(size = 14),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            strip.text = element_text(size = 14, face = "bold"),
            legend.title = element_blank(),
            legend.text = element_text(size = 14),)
    if (tp %in% c("Day29", "Day57")) {
      boxplot_list[[aa]] <- boxplot_list[[aa]] + 
        geom_hline(yintercept = log10(lloqs[assays][aa]), linetype = 2, color = "black", lwd = 1) +
        geom_hline(yintercept = log10(uloqs[assays][aa]), linetype = 2, color = "black", lwd = 1)
    }
    
    if (!is.na(pos.cutoffs[aa])) {
      boxplot_list[[aa]] <- boxplot_list[[aa]] + 
        geom_hline(yintercept = log10(pos.cutoffs[assays][aa]), linetype = 2, color = "black", lwd = 1) 
    } else {
      boxplot_list[[aa]] <- boxplot_list[[aa]] + 
        geom_hline(yintercept = log10(llods[assays][aa]), linetype = 2, color = "black", lwd = 1)
    }
    
    
  }
  
  # Suppress hline warnings
  suppressWarnings(ggsave(ggarrange(plotlist = boxplot_list, 
                   ncol = 2, nrow = 2, 
                   common.legend = TRUE, legend = "bottom",
                   align = "h") + 
           theme(plot.title = element_text(hjust = 0.5, size = 10)),
         filename = paste0(save.results.to, "/boxplots_", tp, "_trt_vaccine_x_cc_",
                           study_name, ".png"),
         height = 9, width = 8))
  
}










##===============================================================================================================
## Box plots, overlayed with boxplots and jittered sample points, for the distribution of baseline, day 57, or 
## fold-rise in assay readouts versus the case-control status among the vaccine recipients, 
## stratified by the baseline serostatus, age group and risk group
##===============================================================================================================
set.seed(12345)
dat.sample5 <-  dat.long.cor.subset %>%
  filter(., .$Trt == "Vaccine") %>%
  split(., list(.$assay, .$cohort_event, .$demo_lab)) %>%
  lapply(., function(x) {
    if(nrow(x) <= 100) {
      return(x)
    } else {
      return(x[sample(1:nrow(x), size = 100),])
    }}) %>% bind_rows



for (tp in tps[tps %in% times]){
  subdat <-dat.long.cor.subset
  subdat_jitter <- dat.sample5
  if (tp %in% c("Day57", "Delta57overB")) {
    subdat <- subdat %>% filter(cohort_event != "Intercurrent Cases" & ph2.D57==1)
    subdat_jitter <- subdat_jitter %>% filter(cohort_event != "Intercurrent Cases" & ph2.D57==1)
  }
  for (aa in seq_along(assays)) {
    boxplots <-  ggplot(subset(subdat, assay == assays[aa] & Trt == "Vaccine"), aes_string(x = "cohort_event", y = tp)) +
      geom_boxplot(aes(colour = cohort_event), width = 0.6, lwd = 1) + 
      stat_boxplot(geom = "errorbar", aes(colour = cohort_event), width = 0.45, lwd = 1) +
      geom_jitter(data = filter(subdat_jitter, assay == assays[aa] & Trt == "Vaccine"),
                  mapping = aes(colour = cohort_event), width = 0.1, 
                  size = 1.4, show.legend = FALSE) +
      theme_pubr() + 
      guides(alpha = "none", color = guide_legend(nrow = 1, byrow = TRUE)) +
      facet_wrap(~ demo_lab) +
      ylab(labels.axis[tp, aa]) + xlab("") + ggtitle(labels.title2[tp, aa]) +
      scale_color_manual(values = c("#1749FF", "#D92321", "#0AB7C9", "#FF6F1B", "#810094", "#378252", 
                                    "#FF5EBF", "#3700A5", "#8F8F8F", "#787873")) +
      scale_y_continuous(limits = c(-2, 6), labels = label_math(10^.x), breaks = seq(-2, 6, 2)) +
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            panel.border = element_rect(fill = NA),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_line(),
            axis.title = element_text(size = 14),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            strip.text = element_text(size = 14, face = "bold"),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 14))
    
    if (tp %in% c("Day29", "Day57")) {
      boxplots <- boxplots + 
        geom_hline(yintercept = log10(llods[assays][aa]), linetype = 2, color = "black", lwd = 1) +
        geom_hline(yintercept = log10(lloqs[assays][aa]), linetype = 2, color = "black", lwd = 1) +
        geom_hline(yintercept = log10(uloqs[assays][aa]), linetype = 2, color = "black", lwd = 1)
    }
    
    
    ggsave(boxplots,
           filename = paste0(save.results.to, "/boxplots_", tp,"_trt_vaccine_x_cc_",
                             assays[aa], "_", study_name, ".png"),
           height = 9,
           width = 8)
  }
}


#-----------------------------------------------
# Spaghetti PLOTS
#-----------------------------------------------
# - Spaghetti plots of antibody marker change over time
#-----------------------------------------------

## in each baseline serostatus group, randomly select 10 placebo recipients and 20 vaccine recipients
set.seed(12345)

# First we need to restrict to ph2.D57==1 if has57==true, to remove participants without day57 ab readings
if(has57) {
  dat.cor.subset.spaghetti <- filter(dat.cor.subset, ph2.D57==1)
  dat.long.cor.subset.spaghetti <- filter(dat.long.cor.subset, ph2.D57==1)
} else {
  dat.cor.subset.spaghetti <- dat.cor.subset
  dat.long.cor.subset.spaghetti <- dat.long.cor.subset
}


var_names <- expand.grid(times = intersect(c("B", "Day29", "Day57"), times),
                         assays = assays) %>%
  mutate(var_names = paste0(times, assays)) %>%
  .[, "var_names"]

spaghetti_ptid <- dat.cor.subset.spaghetti[, c("Ptid", "Trt", var_names, "cohort_event")] %>%
  filter(., complete.cases(.), Trt == 1) %>%
  transmute(cohort_event = cohort_event,
            Ptid = Ptid) %>%
  split(., list(.$cohort_event)) %>%
  lapply(function(xx) {
    if (nrow(xx) <= 20) {
      return(xx$Ptid)
    } else {
      return(sample(xx$Ptid, 20))
    }
  }) %>% unlist %>% as.character

spaghetti_dat <- dat.long.cor.subset.spaghetti[, c("Ptid", "cohort_event", "assay",
                                         intersect(c("B", "Day29", "Day57"), times))] %>%
  filter(Ptid %in% spaghetti_ptid) %>%
  pivot_longer(cols = intersect(c("B", "Day29", "Day57"), times),
               names_to = "time") %>%
  mutate(assay_label = factor(assay, levels = assays, labels = labels.assays.short[assays]),
         time_label = factor(time, levels = intersect(c("B", "Day29", "Day57"), times),
                             labels = c("D1", "D29", "D57")[c("B", "Day29", "Day57") %in% times])) %>%
  as.data.frame


subdat <- spaghetti_dat
covid_corr_spaghetti_facets(plot_dat = subdat,
                            x = "time_label",
                            y = "value",
                            id = "Ptid",
                            color = "cohort_event",
                            facet_by = "assay_label",
                            plot_title = "Vaccine group",
                            ylim = c(-2, 6),
                            ybreaks = seq(-2, 6, by = 2),
                            filename = paste0(
                              save.results.to, "/spaghetti_plot_trt_",
                              study_name, ".png"
                            ),
                            height = 6, width = 5)

