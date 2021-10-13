#Sys.setenv(TRIAL = "moderna_mock") # janssen_pooled_real   janssen_pooled_mock   moderna_mock
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows")
  .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))
#-----------------------------------------------

source(here::here("code", "params.R"))
source(here::here("code", "cor_wrcdf_plot_function.R"))
library(ggpubr)
library(scales)
library(survival)
## get a crude estimate of VE and upper/lower bound of 95% CI
## using an unadjusted Cox model

dat.long.cor.subset <- readRDS(here("data_clean",
                                    "long_cor_data.rds"))

dat.cor.subset <- readRDS(here("data_clean",
                               "cor_data.rds"))

# path for figures
save.results.to <- here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to <- paste0(here::here("output"), "/", attr(config,"config"));
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to <- paste0(save.results.to, "/", COR,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

# Which plots to output
plot_ve_curves <- c(1, 1)
for (bstatus in 0:1) {
    if (verbose) print(bstatus)
  
  # if (has57) {                    ----- NOT NEEDED FOR GENERALIZATION, COMMENTED OUT -----
  #   ve_dat <- dat.cor.subset %>%
  #     filter(EventTimePrimaryD57 >= 7, Bserostatus == bstatus) %>%
  #     transmute(Trt = Trt,
  #               event = EventIndPrimaryD57,
  #               time = EventTimePrimaryD57)
  # 
  #   if (all(ve_dat$event == 0)) {
  #     ## no events in the strata
  #     plot_ve_curves[bstatus + 1] <- 0
  #     print(paste0(
  #       "No D57 events among subjects who are Baseline ",
  #       c("Seronegative", "Seropositive")[bstatus + 1]
  #     ))
  #   } else {
  #     #      if(!file.exists("../cor_coxph/output/D57/marginalized.risk.no.marker."%.%study_name%.%".Rdata")) {
  #     #          # compute ve. this works IF ve_dat is properly defined for computing ve
  #     #          cox_results <- coxph(Surv(time, event) ~ Trt, data = ve_dat)
  #     #          VE <- as.numeric(1 - exp(cox_results$coefficients))
  #     #          VE_lb <- as.numeric(1 - exp(confint.default(cox_results)[2]))
  #     #          VE_ub <- as.numeric(1 - exp(confint.default(cox_results)[1]))
  #     #       } else {
  # 
  #     # use ve computed as part of cor_coxph, which is defined based on marginalized risks
  #     load(paste0(here::here("..", "cor_coxph", "output"), "/",
  #                 attr(config,"config"),
  #                 "/"%.%COR%.%"/marginalized.risk.no.marker."%.%study_name%.%".Rdata"))
  #     VE <- overall.ve[1]
  #     VE_lb <- overall.ve[2]
  #     VE_ub <- overall.ve[3]
  # 
  #     for (aa in 1:length(assays)) {
  #     if (verbose) print(aa)
  #       subdat <-
  #         subset(
  #           dat.long.cor.subset,
  #           assay == assays[aa] &
  #             Trt == "Vaccine" & Bserostatus == bstatus.labels[bstatus + 1]
  #         )
  #       covid_corr_rcdf_ve_lines(
  #         x = subdat$Day57,
  #         weights = subdat$wt.D57,
  #         VE = VE,
  #         VE_lb = VE_lb,
  #         VE_ub = VE_ub,
  #         xlab = labels.axis["Day57", aa],
  #         filename = paste0(
  #           save.results.to,
  #           "/RCDF_VE_Day57_trt_Vaccine_",
  #           bstatus.labels.2[bstatus + 1],
  #           "_",
  #           assays[aa],
  #           "_",
  #           study_name,
  #           ".png"
  #         )
  #       )
  #     }
  # 
  #     # use ve computed as part of cor_coxph, which is defined based on marginalized risks
  #     load(paste0(here::here("..", "cor_coxph", "output"), "/",
  #                 attr(config, "config"),
  #                 "/D29/marginalized.risk.no.marker."%.%study_name%.%".Rdata"))
  #     VE <- overall.ve[1]
  #     VE_lb <- overall.ve[2]
  #     VE_ub <- overall.ve[3]
  # 
  #     for (aa in 1:length(assays)) {
  #       subdat <-
  #         subset(
  #           dat.long.cor.subset,
  #           assay == assays[aa] &
  #             Trt == "Vaccine" & Bserostatus == bstatus.labels[bstatus + 1]
  #         )
  #       covid_corr_rcdf_ve_lines(
  #         x = subdat$Day29,
  #         weights = subdat$wt.D29,
  #         VE = VE,
  #         VE_lb = VE_lb,
  #         VE_ub = VE_ub,
  #         xlab = labels.axis["Day29", aa],
  #         filename = paste0(
  #           save.results.to,
  #           "/RCDF_VE_Day29_trt_Vaccine_",
  #           bstatus.labels.2[bstatus + 1],
  #           "_",
  #           assays[aa],
  #           "_",
  #           study_name,
  #           ".png"
  #         )
  #       )
  #     }
  #   }
  # } else {
  
    ve_dat <- dat.cor.subset %>%
      filter(!!as.symbol(config.cor$EventTimePrimary) >= config.cor$tpeaklag, Bserostatus == bstatus) %>%
      transmute(Trt = Trt,
                event = !!as.symbol(config.cor$EventIndPrimary),
                time = !!as.symbol(config.cor$EventTimePrimary))

    if (all(ve_dat$event == 0)) {
      ## no events in the strata
      plot_ve_curves[bstatus + 1] <- 0
      print(paste0(
        "No "%.%COR%.%" events among subjects who are Baseline ",
        c("Seronegative", "Seropositive")[bstatus + 1]
      ))
    } else {
      
      if(!file.exists(paste0(here::here("..", "cor_coxph", "output"), "/", attr(config, "config"), "/"%.%COR%.%"/marginalized.risk.no.marker."%.%study_name%.%".Rdata"))) {
         # compute ve. this works IF ve_dat is properly defined for computing ve
         cox_results <- coxph(Surv(time, event) ~ Trt, data = ve_dat)
         VE <- as.numeric(1 - exp(cox_results$coefficients))
         VE_lb <- as.numeric(1 - exp(confint.default(cox_results)[2]))
         VE_ub <- as.numeric(1 - exp(confint.default(cox_results)[1]))
         
      } else {
      
      # use ve computed as part of cor_coxph, which is defined based on marginalized risks
      load(paste0(here::here("..", "cor_coxph", "output"), "/",
                  attr(config, "config"),
                  "/"%.%COR%.%"/marginalized.risk.no.marker."%.%study_name%.%".Rdata"))
      VE <- overall.ve[1]
      VE_lb <- overall.ve[2]
      VE_ub <- overall.ve[3]
      
      }

      for (aa in 1:length(assays)) {
        subdat <-
          subset(
            dat.long.cor.subset,
            assay == assays[aa] &
              Trt == "Vaccine" & Bserostatus == bstatus.labels[bstatus + 1]
          )
        covid_corr_rcdf_ve_lines(
          x = subdat[[paste0("Day", tpeak)]],
          weights = subdat[[config.cor$wt]],
          VE = VE,
          VE_lb = VE_lb,
          VE_ub = VE_ub,
          xlab = labels.axis[paste0("Day", tpeak), aa],
          filename = paste0(
            save.results.to,
            "/RCDF_VE_", "Day", tpeak, "_trt_Vaccine_",
            bstatus.labels.2[bstatus + 1],
            "_",
            assays[aa],
            "_",
            study_name,
            ".png"
          )
        )
      }
    }
  # }
}

saveRDS(plot_ve_curves,
        here::here("data_clean", "plot_ve_curves.rds"))
