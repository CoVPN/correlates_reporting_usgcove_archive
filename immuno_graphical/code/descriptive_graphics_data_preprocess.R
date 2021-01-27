#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(dplyr)
library(here)
library(stringr)
library(COVIDcorr)

# load data
data(dat.mock)

# load parameters 
source(here("code", "params.R"))

## setting the floor values
dat <- dat.mock %>%
  mutate(
    BbindSpike = ifelse(BbindSpike >= log10(17), BbindSpike, log10(17)),
    Day29bindSpike = ifelse(Day29bindSpike >= log10(17), Day29bindSpike,
                            log10(17)),
    Day57bindSpike = ifelse(Day57bindSpike >= log10(17), Day57bindSpike,
                            log10(17)),
    BbindRBD = ifelse(BbindRBD >= log10(17), BbindRBD, log10(17)),
    Day29bindRBD = ifelse(Day29bindRBD >= log10(17), Day29bindRBD, log10(17)),
    Day57bindRBD = ifelse(Day57bindRBD >= log10(17), Day57bindRBD, log10(17)),
    Bpseudoneutid50 = ifelse(Bpseudoneutid50 >= log10(25), Bpseudoneutid50,
                             log10(25)),
    Day29pseudoneutid50 = ifelse(Day29pseudoneutid50 >= log10(25),
                                 Day29pseudoneutid50, log10(25)),
    Day57pseudoneutid50 = ifelse(Day57pseudoneutid50 >= log10(25),
                                 Day57pseudoneutid50, log10(25)),
    Bpseudoneutid80 = ifelse(Bpseudoneutid80 >= log10(22), Bpseudoneutid80,
                             log10(22)),
    Day29pseudoneutid80 = ifelse(Day29pseudoneutid80 >= log10(22),
                                 Day29pseudoneutid80, log10(22)),
    Day57pseudoneutid80 = ifelse(Day57pseudoneutid80 >= log10(22),
                                 Day57pseudoneutid80, log10(22))
  )


# For immunogenicity characterization, complete ignore any information on cases
# vs. non-cases.  The goal is to characterize immunogenicity in the random
# subcohort, which is a stratified sample of enrolled participants. So,
# immunogenicity analysis is always done in ppts that meet all of the criteria.
dat.twophase.sample <- dat %>%
  dplyr::filter(TwophasesampInd== 1 & SubcohortInd == 1 & Perprotocol == 1)
twophase_sample_id <- dat.twophase.sample$Ptid

# First focus on baseline negative pooling over all baseline demog strata. This
# is the primary cohort for CoR analysis. So, all plots and tables relative to
# this could be shown first in the pdf. Second output contrasting results
# baseline negative vs. baseline positive, again pooling over all baseline
# demog strata.
# Supp material that expands
#  1. for the individual baseline demog cells, for completeness.
#  2. for the individual baseline demog cells, for completeness.
dat$EventLabelD29 <- factor(dat$EventIndPrimaryD29, levels = c(0, 1),
                            labels = c("D29 Non-Case", "D29 Case"))
dat$EventLabelD57 <- factor(dat$EventIndPrimaryD57, levels = c(0, 1),
                            labels = c("D57 Non-Case", "D57 Case"))

## arrange the dataset in the long form, expand by assay types
## dat.long.1 is the subject level covariates;
## dat.long.2 is the long-form time varying variables
dat.long.1.0 <- dat[, c("Ptid", "Trt", "MinorityInd", "HighRiskInd", "Age",
                        "Sex", "Bserostatus", "Fullvaccine", "Perprotocol",
                        "EventIndPrimaryD29", "EventIndPrimaryD57",
                        "SubcohortInd", "age.geq.65",  "TwophasesampInd",
                        "Bstratum", "wt", "EventLabelD29", "EventLabelD57",
                        "race", "ethnicity", "WhiteNonHispanic")]

dat.long.1 <- bind_rows(replicate(4, dat.long.1.0, simplify = FALSE))
name_grid <- expand.grid(aa = times,
                         cc = c("", "CPV", paste(".imp", 1:10, sep = "")))
dat.long.2.names <- paste(name_grid$aa, name_grid$cc, sep = "")
dat.long.2 <- as.data.frame(matrix(nrow = nrow(dat) * 4,
                                   ncol = length(dat.long.2.names)))
colnames(dat.long.2) <- dat.long.2.names

for (ii in 1:nrow(name_grid)) {
  dat_mock_col_names <- paste(name_grid$aa[ii], assays,
                              name_grid$cc[ii], sep = "")
  dat.long.2[, dat.long.2.names[ii]] <-
    unlist(lapply(dat_mock_col_names, function(nn) {
                    if (nn %in% colnames(dat)) {
                      dat[, nn]
                    } else {
                      rep(NA, nrow(dat))
                    }
           })
    )
}

dat.long.2$assay <- rep(assays, each = nrow(dat))
dat.long <- cbind(dat.long.1, dat.long.2)

## change the labels of the factors for plot labels
dat.long$Trt <- factor(dat.long$Trt, levels = c(0, 1), labels = trt.labels)
dat.long$Bserostatus <- factor(dat.long$Bserostatus, levels = c(0, 1),
                               labels = bstatus.labels)
dat.long$assay <- factor(dat.long$assay, levels = assays, labels = assays)

dat.long.twophase.sample <- dat.long[dat.long$Ptid %in% twophase_sample_id, ]
dat.twophase.sample <- subset(dat, Ptid %in% twophase_sample_id)

## label the subjects according to their case-control status
dat.long.twophase.sample <- dat.long.twophase.sample %>%
  mutate(
    EventD29 = factor(EventIndPrimaryD29, levels = c(0, 1),
                      labels = c("Non-Case", "Case")),
    EventD57 = factor(EventIndPrimaryD57, levels = c(0, 1),
                      labels = c("Non-Case", "Case"))
  )
dat.long <- dat.long %>%
  mutate(
    EventD29 = factor(EventIndPrimaryD29, levels = c(0, 1),
                      labels = c("Non-Case", "Case")),
    EventD57 = factor(EventIndPrimaryD57, levels = c(0, 1),
                      labels = c("Non-Case", "Case"))
  )

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
        file = here("data_clean", "long_twophase_data.rds"))
saveRDS(dat.twophase.sample,
        file = here("data_clean", "twophase_data.rds"))
