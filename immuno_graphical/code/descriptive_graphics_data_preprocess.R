#--------------------------------------------
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#--------------------------------------------

library(dplyr)
library(here)
library(stringr)
library(COVIDcorr)

dat <- COVIDcorr::dat.mock$dat.mock

source(here::here("code", "params.R"))

## setting the floor values
dat$BbindSpike <- ifelse(dat$BbindSpike >= log10(17), dat$BbindSpike, log10(17))
dat$Day29bindSpike <- ifelse(dat$Day29bindSpike >= log10(17), dat$Day29bindSpike, log10(17))
dat$Day57bindSpike <- ifelse(dat$Day57bindSpike >= log10(17), dat$Day57bindSpike, log10(17))

dat$BbindRBD <- ifelse(dat$BbindRBD >= log10(17), dat$BbindRBD, log10(17))
dat$Day29bindRBD <- ifelse(dat$Day29bindRBD >= log10(17), dat$Day29bindRBD, log10(17))
dat$Day57bindRBD <- ifelse(dat$Day57bindRBD >= log10(17), dat$Day57bindRBD, log10(17))

dat$Bpseudoneutid50 <- ifelse(dat$Bpseudoneutid50 >= log10(25), dat$Bpseudoneutid50, log10(25))
dat$Day29pseudoneutid50 <- ifelse(dat$Day29pseudoneutid50 >= log10(25), dat$Day29pseudoneutid50, log10(25))
dat$Day57pseudoneutid50 <- ifelse(dat$Day57pseudoneutid50 >= log10(25), dat$Day57pseudoneutid50, log10(25))

dat$Bpseudoneutid80 <- ifelse(dat$Bpseudoneutid80 >= log10(22), dat$Bpseudoneutid80, log10(22))
dat$Day29pseudoneutid80 <- ifelse(dat$Day29pseudoneutid80 >= log10(22), dat$Day29pseudoneutid80, log10(22))
dat$Day57pseudoneutid80 <- ifelse(dat$Day57pseudoneutid80 >= log10(22), dat$Day57pseudoneutid80, log10(22))

## add the floor value of Live Virus Neutralization Titer when the data are available
# dat$Bliveneutid50 <- ifelse(dat$Bliveneutid50 >= log10(25), dat$Bliveneutid50, log10(25))
# dat$Day29liveneutid50 <- ifelse(dat$Day29liveneutid50 >= log10(25), dat$Day29liveneutid50, log10(25))
# dat$Day57liveneutid50 <- ifelse(dat$Day57liveneutid50 >= log10(25), dat$Day57liveneutid50, log10(25))
# 
# dat$Bliveneutid80 <- ifelse(dat$Bliveneutid80 >= log10(22), dat$Bliveneutid80, log10(22))
# dat$Day29liveneutid80 <- ifelse(dat$Day29liveneutid80 >= log10(22), dat$Day29liveneutid80, log10(22))
# dat$Day57liveneutid80 <- ifelse(dat$Day57liveneutid80 >= log10(22), dat$Day57liveneutid80, log10(22))


## For immunogenicity characterization, complete ignore any information on cases vs. non-cases.  The goal is to
## characterize immunogenicity in the random subcohort, which is a stratified sample of enrolled participants.
## So immunogenicity analysis is always done in ppts that meet all of the criteria

dat.twophase.sample <- filter(dat, TwophasesampInd== 1 & SubcohortInd == 1 & Perprotocol == 1)


twophase_sample_id <- dat.twophase.sample$Ptid

#First focus on baseline negative pooling over all baseline demog strata.  This is the primary cohort for CoR analysis.  So all plots / tables relative to this could be shown first in the pdf.
#Second output contrasting results baseline negative vs. baseline positive, again pooling over all baseline demog strata.
#Supp material that expands 1. for the individual baseline demog cells, for completeness.
#Supp material that expands 2. for the individual baseline demog cells, for completeness.

dat$EventLabelD29 <- factor(dat$EventIndPrimaryD29, levels = c(0, 1), labels = c("D29 Non-Case", "D29 Case"))
dat$EventLabelD57 <- factor(dat$EventIndPrimaryD57, levels = c(0, 1), labels = c("D57 Non-Case", "D57 Case"))

## arrange the dataset in the long form, expand by assay types
## dat.long.subject_level is the subject level covariates;
## dat.long.assay_value is the long-form time variables that differ by the assay type
dat.long.subject_level <- dat[, c("Ptid", "Trt", "MinorityInd", "HighRiskInd", "Age", "Sex",
                        "Bserostatus", "Fullvaccine", "Perprotocol", "EventIndPrimaryD29",
                        "EventIndPrimaryD57", "SubcohortInd", "age.geq.65",  "TwophasesampInd", 
                        "Bstratum", "wt", "EventLabelD29", "EventLabelD57", "race", "ethnicity", 
                        "WhiteNonHispanic")] %>%
  replicate(length(assays), ., simplify = FALSE) %>%
  bind_rows


name_grid <- expand.grid(aa = times,
                         cc = c("", "CPV", paste(".imp", 1:10, sep = "")))
dat.long.assay_value.names <- paste(name_grid$aa, name_grid$cc, sep = "")
dat.long.assay_value <- as.data.frame(matrix(nrow = nrow(dat) * length(assays),
                                        ncol = length(dat.long.assay_value.names)))
colnames(dat.long.assay_value) <- dat.long.assay_value.names

for (ii in 1:nrow(name_grid)) {
  dat_mock_col_names <- paste(name_grid$aa[ii], assays, name_grid$cc[ii], sep = "")
  dat.long.assay_value[, dat.long.assay_value.names[ii]] <- unlist(lapply(dat_mock_col_names,
                                                                function(nn) {
                                                                  if (nn %in% colnames(dat)) {
                                                                    dat[, nn]
                                                                  } else {
                                                                    rep(NA, nrow(dat))
                                                                  }
                                                                }))

}

dat.long.assay_value$assay <- rep(assays, each = nrow(dat))

dat.long <- cbind(dat.long.subject_level, dat.long.assay_value)

## change the labels of the factors for plot labels
dat.long$Trt <- factor(dat.long$Trt, levels = c(0, 1), labels = trt.labels)
dat.long$Bserostatus <- factor(dat.long$Bserostatus, levels = c(0, 1), labels = bstatus.labels)
dat.long$assay <- factor(dat.long$assay, levels = assays,
                              labels = assays)

dat.long.twophase.sample <- dat.long[dat.long$Ptid %in% twophase_sample_id, ]
dat.twophase.sample <- subset(dat, Ptid %in% twophase_sample_id)

## label the subjects according to their case-control status
dat.long.twophase.sample$EventD29 <- factor(dat.long.twophase.sample$EventIndPrimaryD29, levels = c(0, 1), labels = c("Non-Case", "Case"))
dat.long.twophase.sample$EventD57 <- factor(dat.long.twophase.sample$EventIndPrimaryD57, levels = c(0, 1), labels = c("Non-Case", "Case"))
dat.long$EventD29 <- factor(dat.long$EventIndPrimaryD29, levels = c(0, 1), labels = c("Non-Case", "Case"))
dat.long$EventD57 <- factor(dat.long$EventIndPrimaryD57, levels = c(0, 1), labels = c("Non-Case", "Case"))



# # matrix to decide the sampling strata
dat.long$demo_lab <-
  with(dat.long, factor(paste0(age.geq.65, HighRiskInd),
                        levels = c("00", "01", "10", "11"),
                        labels = c("Age < 65 not at tisk",
                                   "Age < 65 at risk",
                                   "Age >= 65 not at risk",
                                   "Age >= 65 at risk")))

# labels of the demographic strata for the subgroup plotting
dat.long.twophase.sample$trt_bstatus_label <- 
  with(dat.long.twophase.sample, 
       factor(paste0(as.numeric(Trt), as.numeric(Bserostatus)),
              levels = c("11", "12", "21", "22"),
              labels = c("Placebo, Baseline Neg",
                         "Placebo, Baseline Pos",
                         "Vaccine, Baseline Neg",
                         "Vaccine, Baseline Pos")))



dat.long.twophase.sample$age_geq_65_label <- 
  with(dat.long.twophase.sample,
       factor(age.geq.65, 
              levels = c(0, 1), 
              labels = c("Age >= 65", "Age < 65")))

dat.long.twophase.sample$highrisk_label <- 
  with(dat.long.twophase.sample,
       factor(HighRiskInd, 
              levels = c(0, 1), 
              labels = c("Not at risk", "At risk")))

dat.long.twophase.sample$age_risk_label <- 
  with(dat.long.twophase.sample, 
       factor(paste0(age.geq.65, HighRiskInd),
              levels = c("00", "01", "10", "11"),
              labels = c("Age < 65 not at tisk",
                         "Age < 65 at risk",
                         "Age >= 65 not at risk",
                         "Age >= 65 at risk")))

dat.long.twophase.sample$sex_label <- 
  with(dat.long.twophase.sample,
       factor(Sex,
              levels = c(0, 1), 
              labels = c("Female", "Male")))

dat.long.twophase.sample$age_sex_label <- 
  with(dat.long.twophase.sample, 
       factor(paste0(age.geq.65, Sex),
              levels = c("00", "01", "10", "11"),
              labels = c("Age < 65 male",
                         "Age < 65 female",
                         "Age >= 65 male",
                         "Age >= 65 female")))

dat.long.twophase.sample$minority_label <-
  with(dat.long.twophase.sample,
       factor(MinorityInd, 
              levels = c(0, 1),
              labels = c("White Non-Hispanic", "Comm. of Color")))

dat.long.twophase.sample$age_minority_label <- 
  with(dat.long.twophase.sample, 
       factor(paste0(age.geq.65, MinorityInd),
              levels = c("00", "01", "10", "11"),
              labels = c("Age < 65 Comm. of Color",
                         "Age < 65 White Non-Hispanic",
                         "Age >= 65 Comm. of Color",
                         "Age >= 65 White Non-Hispanic")))

saveRDS(dat.long.twophase.sample, 
        file = here::here("data_clean", "long_twophase_data.rds"))
saveRDS(dat.twophase.sample, 
        file = here::here("data_clean", "twophase_data.rds"))
