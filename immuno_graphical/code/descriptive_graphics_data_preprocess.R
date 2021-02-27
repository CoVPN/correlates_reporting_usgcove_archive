#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

### reset the floor values
### set the ULOQ values
### set the ethnicity value
### add 



library(here)
library(dplyr)
library(stringr)
# DB: Scheduled for deletion
# library(COVIDcorr)
# # load data
# data(dat.mock)

dat.mock <- read.csv(here("..", "data_raw", data_name), header = TRUE)

# load parameters
source(here("code", "params.R"))


dat <- dat.mock 

## setting the floor values
dat <- dat %>% mutate(
    BbindSpike = ifelse(BbindSpike >= log10(20), BbindSpike, log10(10)),
    Day29bindSpike = ifelse(Day29bindSpike >= log10(20), Day29bindSpike,
      log10(10)
    ),
    Day57bindSpike = ifelse(Day57bindSpike >= log10(20), Day57bindSpike,
      log10(10)
    ),
    BbindRBD = ifelse(BbindRBD >= log10(20), BbindRBD, log10(10)),
    Day29bindRBD = ifelse(Day29bindRBD >= log10(20), Day29bindRBD, log10(10)),
    Day57bindRBD = ifelse(Day57bindRBD >= log10(20), Day57bindRBD, log10(10)),
    Bpseudoneutid50 = ifelse(Bpseudoneutid50 >= log10(20), Bpseudoneutid50,
      log10(10)
    ),
    Day29pseudoneutid50 = ifelse(Day29pseudoneutid50 >= log10(20),
      Day29pseudoneutid50, log10(10)
    ),
    Day57pseudoneutid50 = ifelse(Day57pseudoneutid50 >= log10(20),
      Day57pseudoneutid50, log10(10)
    ),
    Bpseudoneutid80 = ifelse(Bpseudoneutid80 >= log10(20), Bpseudoneutid80,
      log10(10)
    ),
    Day29pseudoneutid80 = ifelse(Day29pseudoneutid80 >= log10(20),
      Day29pseudoneutid80, log10(10)
    ),
    Day57pseudoneutid80 = ifelse(Day57pseudoneutid80 >= log10(20),
      Day57pseudoneutid80, log10(10)
    )
  )


# For immunogenicity characterization, complete ignore any information on cases
# vs. non-cases.  The goal is to characterize immunogenicity in the random
# subcohort, which is a stratified sample of enrolled participants. So,
# immunogenicity analysis is always done in ppts that meet all of the criteria.
dat.twophase.sample <- dat %>%
  dplyr::filter(TwophasesampInd == 1 & SubcohortInd == 1 & Perprotocol == 1)
twophase_sample_id <- dat.twophase.sample$Ptid


## arrange the dataset in the long form, expand by assay types
## dat.long.subject_level is the subject level covariates;
## dat.long.assay_value is the long-form time variables that differ by the assay type
dat.long.subject_level <- dat[, c(
  "Ptid", "Trt", "MinorityInd", "HighRiskInd", "Age", "Sex",
  "Bserostatus", "Fullvaccine", "Perprotocol", "EventIndPrimaryD29",
  "EventIndPrimaryD57", "SubcohortInd", "age.geq.65", "TwophasesampInd",
  "Bstratum", "wt", "race",
  "EthnicityHispanic","EthnicityNotreported", "EthnicityUnknown"
  "WhiteNonHispanic"
)] %>%
  replicate(length(assays), ., simplify = FALSE) %>%
  bind_rows()

name_grid <- expand.grid(
  aa = times,
  cc = c("", "CPV", paste(".imp", 1:10, sep = ""))
)

dat.long.assay_value.names <- paste(name_grid$aa, name_grid$cc, sep = "")
dat.long.assay_value <- as.data.frame(matrix(
  nrow = nrow(dat) * length(assays),
  ncol = length(dat.long.assay_value.names)
))
colnames(dat.long.assay_value) <- dat.long.assay_value.names

for (ii in 1:nrow(name_grid)) {
  dat_mock_col_names <- paste(name_grid$aa[ii], assays, name_grid$cc[ii], sep = "")
  dat.long.assay_value[, dat.long.assay_value.names[ii]] <- unlist(lapply(
    dat_mock_col_names,
    function(nn) {
      if (nn %in% colnames(dat)) {
        dat[, nn]
      } else {
        rep(NA, nrow(dat))
      }
    }
  ))
}

dat.long.assay_value$assay <- rep(assays, each = nrow(dat))

dat.long <- cbind(dat.long.subject_level, dat.long.assay_value)


## change the labels of the factors for plot labels
dat.long$Trt <- factor(dat.long$Trt, levels = c(0, 1), labels = trt.labels)
dat.long$Bserostatus <- factor(dat.long$Bserostatus,
  levels = c(0, 1),
  labels = bstatus.labels
)
dat.long$assay <- factor(dat.long$assay, levels = assays, labels = assays)

dat.long.twophase.sample <- dat.long[dat.long$Ptid %in% twophase_sample_id, ]
dat.twophase.sample <- subset(dat, Ptid %in% twophase_sample_id)

## label the subjects according to their case-control status
dat.long.twophase.sample <- dat.long.twophase.sample %>%
  mutate(
    EventD29 = factor(EventIndPrimaryD29,
      levels = c(0, 1),
      labels = c("Non-Case", "Case")
    ),
    EventD57 = factor(EventIndPrimaryD57,
      levels = c(0, 1),
      labels = c("Non-Case", "Case")
    )
  )
dat.long <- dat.long %>%
  mutate(
    EventD29 = factor(EventIndPrimaryD29,
      levels = c(0, 1),
      labels = c("Non-Case", "Case")
    ),
    EventD57 = factor(EventIndPrimaryD57,
      levels = c(0, 1),
      labels = c("Non-Case", "Case")
    )
  )

# # matrix to decide the sampling strata
dat.long$demo_lab <-
  with(dat.long, factor(paste0(age.geq.65, HighRiskInd),
    levels = c("00", "01", "10", "11"),
    labels = c(
      "Age < 65 not at risk",
      "Age < 65 at risk",
      "Age >= 65 not at risk",
      "Age >= 65 at risk"
    )
  ))

# labels of the demographic strata for the subgroup plotting
dat.long.twophase.sample$trt_bstatus_label <-
  with(
    dat.long.twophase.sample,
    factor(paste0(as.numeric(Trt), as.numeric(Bserostatus)),
      levels = c("11", "12", "21", "22"),
      labels = c(
        "Placebo, Baseline Neg",
        "Placebo, Baseline Pos",
        "Vaccine, Baseline Neg",
        "Vaccine, Baseline Pos"
      )
    )
  )

dat.long.twophase.sample$age_geq_65_label <-
  with(
    dat.long.twophase.sample,
    factor(age.geq.65,
      levels = c(0, 1),
      labels = c("Age >= 65", "Age < 65")
    )
  )

dat.long.twophase.sample$highrisk_label <-
  with(
    dat.long.twophase.sample,
    factor(HighRiskInd,
      levels = c(0, 1),
      labels = c("Not at risk", "At risk")
    )
  )

dat.long.twophase.sample$age_risk_label <-
  with(
    dat.long.twophase.sample,
    factor(paste0(age.geq.65, HighRiskInd),
      levels = c("00", "01", "10", "11"),
      labels = c(
        "Age < 65 not at risk",
        "Age < 65 at risk",
        "Age >= 65 not at risk",
        "Age >= 65 at risk"
      )
    )
  )

dat.long.twophase.sample$sex_label <-
  with(
    dat.long.twophase.sample,
    factor(Sex,
      levels = c(0, 1),
      labels = c("Female", "Male")
    )
  )

dat.long.twophase.sample$age_sex_label <-
  with(
    dat.long.twophase.sample,
    factor(paste0(age.geq.65, Sex),
      levels = c("00", "01", "10", "11"),
      labels = c(
        "Age < 65 male",
        "Age < 65 female",
        "Age >= 65 male",
        "Age >= 65 female"
      )
    )
  )


dat.long.twophase.sample$ethnicity_label <-
  with(
    dat.long.twophase.sample,
    ifelse(EthnicityHispanic == 1,
           "Hispanic or Latino",
           ifelse(
             EthnicityNotreported == 0 & EthnicityUnknown == 0,
             "Not Hispanic or Latino",
             "Not reported and unknown"
           ))
  ) %>% factor(
    levels = c("Hispanic or Latino", "Not Hispanic or Latino", "Not reported and unknown")
  )



dat.long.twophase.sample$minority_label <-
  with(
    dat.long.twophase.sample,
    factor(WhiteNonHispanic,
      levels = c(1, 0),
      labels = c("White Non-Hispanic", "Comm. of Color")
    )
  )

dat.long.twophase.sample$age_minority_label <-
  with(
    dat.long.twophase.sample,
    factor(paste0(age.geq.65, WhiteNonHispanic),
      levels = c("00", "01", "10", "11"),
      labels = c(
        "Age < 65 Comm. of Color",
        "Age < 65 White Non-Hispanic",
        "Age >= 65 Comm. of Color",
        "Age >= 65 White Non-Hispanic"
      )
    )
  )
dat.long.twophase.sample$ethnicity <- as.factor(dat.long.twophase.sample$ethnicity)
dat.twophase.sample$ethnicity <- as.factor(dat.twophase.sample$ethnicity)
dat.long.twophase.sample$race <- as.factor(dat.long.twophase.sample$race)
dat.twophase.sample$race <- as.factor(dat.twophase.sample$race)


saveRDS(as.data.frame(dat.long.twophase.sample),
  file = here("data_clean", "long_twophase_data.rds")
)
saveRDS(as.data.frame(dat.twophase.sample),
  file = here("data_clean", "twophase_data.rds")
)


