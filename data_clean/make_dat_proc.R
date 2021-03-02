#-----------------------------------------------
renv::activate()

# # obligatory to append to the top of each script
# # There is a bug on Windows that prevents renv from working properly. saved.system.libPaths provides a workaround:
# if (.Platform$OS.type == "windows") saved.system.libPaths=.libPaths()
# renv::activate(project = here::here(".."))
# if (.Platform$OS.type == "windows") {
#     options(renv.config.install.transactional = FALSE)
#     renv::restore(library=saved.system.libPaths, prompt=FALSE) # for a quick test, add: packages="backports"
#     .libPaths(c(saved.system.libPaths, .libPaths()))
# } else renv::activate(project = here::here(".."))

# packages and settings
library(here)
library(tidyverse)
library(Hmisc) # wtd.quantile, cut2
library(mice)
library(dplyr)

source(here::here("_common.R"))
#-----------------------------------------------

# load data and rename first column (ID)
dat_proc <- read_csv(here(
  "data_raw", data_in_file
))
colnames(dat_proc)[1] <- "Ptid"

# define age cutoff and two-phase sampling indicator
dat_proc <- dat_proc %>%
  mutate(
    age.geq.65 = as.integer(Age >= 65),
    # create two-phase sampling indicator
    TwophasesampInd = Fullvaccine == 1 &
      (SubcohortInd | EventIndPrimaryD29 == 1) &
      complete.cases(cbind(
        BbindSpike, Day29bindSpike, Day57bindSpike,
        BbindRBD, Day29bindRBD, Day57bindRBD
      ))
  )

# ethnicity labeling
dat_proc$ethnicity <- ifelse(dat_proc$EthnicityHispanic == 1,
  labels.ethnicity[1], labels.ethnicity[2]
)
dat_proc$ethnicity[dat_proc$EthnicityNotreported == 1 |
  dat_proc$EthnicityUnknown == 1] <- labels.ethnicity[3]
dat_proc$ethnicity <- factor(dat_proc$ethnicity, levels = labels.ethnicity)

# race labeling
dat_proc <- dat_proc %>%
  mutate(
    race = labels.race[1],
    race = case_when(
      Black == 1 ~ labels.race[2],
      Asian == 1 ~ labels.race[3],
      NatAmer == 1 ~ labels.race[4],
      PacIsl == 1 ~ labels.race[5],
      Multiracial == 1 ~ labels.race[6],
      Other == 1 ~ labels.race[7],
      Notreported == 1 | Unknown == 1 ~ labels.race[8],
      TRUE ~ labels.race[1]
    ),
    race = factor(race, levels = labels.race)
  )

# WhiteNonHispanic=1 IF race is White AND ethnicity is not Hispanic
dat_proc$WhiteNonHispanic <- NA
dat_proc$WhiteNonHispanic <-
  ifelse(dat_proc$race == "White" &
    dat_proc$ethnicity == "Not Hispanic or Latino", 1,
  dat_proc$WhiteNonHispanic
  )
# WhiteNonHispanic=0 IF race is not "white or unknown" OR ethnicity is Hispanic
dat_proc$WhiteNonHispanic <-
  ifelse(!dat_proc$race %in% c(labels.race[1], last(labels.race)) |
    dat_proc$ethnicity == "Hispanic or Latino", 0,
    dat_proc$WhiteNonHispanic
  )

# check coding via tables
#table(dat_proc$race, useNA = "ifany")
#table(dat_proc$WhiteNonHispanic, useNA = "ifany")
#table(dat_proc$race, dat_proc$WhiteNonHispanic, useNA = "ifany")


###############################################################################
# stratum variables
###############################################################################

# For Moderna
# Bstratum, 1-3, defines the 3 baseline strata within trt/serostatus
dat_proc <- dat_proc %>%
  mutate(
    Bstratum = ifelse(Age >= 65, 1, ifelse(HighRiskInd == 1, 2, 3))
  )
names(Bstratum.labels) <- Bstratum.labels

# tps stratum, 1-12, used in tps regression
dat_proc <- dat_proc %>%
  mutate(
    tps.stratum = Bstratum + strtoi(paste0(Trt, Bserostatus), base = 2) *
      length(Bstratum.labels)
  )

# Wstratum, 1-13. Differs from tps stratum in that case is a separate stratum;
# used to compute sampling weights. NOTE: the case is defined using D29 status
dat_proc <- dat_proc %>%
  mutate(
    Wstratum = ifelse(EventIndPrimaryD29 == 1, 1 + max(tps.stratum),
      tps.stratum
    )
  )


###############################################################################
# observation-level weights
###############################################################################

wts_table <- with(
  subset(dat_proc, Perprotocol == 1 & EventTimePrimaryD57>=7),
  table(Wstratum, TwophasesampInd)
)
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_proc$wt <- wts_norm[dat_proc$Wstratum]
dat_proc$wt[!with(dat_proc, Perprotocol == 1 & EventTimePrimaryD57>=7)] <- NA

wts_table2 <- with(
  subset(dat_proc, EventTimePrimaryD29>=14 & Perprotocol == 1 | EventTimePrimaryD29>=7 & EventTimePrimaryD29<=13 & Fullvaccine==1),
  table(Wstratum, TwophasesampInd)
)
wts_norm2 <- rowSums(wts_table2) / wts_table2[, 2]
dat_proc$wt.2 <- wts_norm2[dat_proc$Wstratum]
dat_proc$wt.2[!with(dat_proc, EventTimePrimaryD29>=14 & Perprotocol == 1 | EventTimePrimaryD29>=7 & EventTimePrimaryD29<=13 & Fullvaccine==1)] <- NA


###############################################################################
# impute missing neut biomarkers in ph2
#     impute vaccine and placebo separately
#     use all assays
#     use baseline and D57, but not Delta
###############################################################################

n.imp <- 10
dat.tmp.impute <- subset(dat_proc, TwophasesampInd == 1)

for (trt in unique(dat_proc$Trt)) {
  imp <- dat.tmp.impute %>%
    dplyr::filter(Trt == trt) %>%
    select(all_of(markers)) %>%
    mice(m = n.imp, printFlag = FALSE)

  dat.tmp.impute[dat.tmp.impute$Trt == trt, markers] <-
    mice::complete(imp, action = 1)
}

stopifnot(
  all(table(dat.tmp.impute$Wstratum, complete.cases(dat.tmp.impute[markers])))
)

# populate dat_proc markers with the imputed values
dat_proc[markers] <-
  dat.tmp.impute[markers][match(dat_proc$Ptid, dat.tmp.impute$Ptid), ]


###############################################################################
# define delta for dat_proc
###############################################################################

dat_proc["Delta57overB" %.% assays] <-
  dat_proc["Day57" %.% assays] - dat_proc["B" %.% assays]
dat_proc["Delta29overB" %.% assays] <-
  dat_proc["Day29" %.% assays] - dat_proc["B" %.% assays]
dat_proc["Delta57over29" %.% assays] <-
  dat_proc["Day57" %.% assays] - dat_proc["Day29" %.% assays]


###############################################################################
# censoring values below LLOD
# llods defined in _common.R
###############################################################################

for (a in assays) {
  for (t in times) {
    dat_proc[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] < log10(llods[a]),
                                  log10(llods[a] / 2), dat_proc[[t %.% a]])
  }
}


###############################################################################
# bundle data sets and save as CSV
###############################################################################
write_csv(dat_proc, file = here("data_clean", data_name))