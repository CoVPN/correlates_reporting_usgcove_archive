#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

source(here::here("code", "cor_data_preprocess.R"))
source(here::here("code", "cor_wrcdf_plot_function.R"))
library(ggpubr)
library(scales)
library(survival)
## get a crude estimate of VE and upper/lower bound of 95% CI 
## using an unadjusted Cox model
for (bstatus in 0:1) {
  ve_dat <- dat.cor.subset %>%
    filter(EventTimePrimaryD57 >= 7, Bserostatus == bstatus) %>%
    transmute(Trt = Trt,
              event = EventIndPrimaryD57,
              time = EventTimePrimaryD57)
  
  cox_results <- coxph(Surv(time, event) ~ Trt, data = ve_dat)
  
  VE <- as.numeric(1 - exp(cox_results$coefficients))
  VE_lb <- as.numeric(1 - exp(confint.default(cox_results)[2]))
  VE_ub <- as.numeric(1 - exp(confint.default(cox_results)[1]))
  
  for (aa in 1:length(assays)) {
    subdat <- subset(dat.long.cor.subset, assay == assays[aa] & Trt == "Vaccine" & Bserostatus == bstatus.labels[bstatus + 1])
    covid_corr_rcdf_ve_lines(x = subdat$Day57,
                             weights = subdat$wt,
                             VE = VE,
                             VE_lb = VE_lb,
                             VE_ub = VE_ub,
                             xlab = labels.axis["Day57", aa],
                             filename = paste0(save.results.to, "/RCDF_VE_Day57_trt_Vaccine_", 
                                               bstatus.labels.2[bstatus + 1], "_", assays[aa], "_", study.name, ".png"))
  }
}
