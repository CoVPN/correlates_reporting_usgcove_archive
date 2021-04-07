#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(dplyr)
library(stringr)
dat.mock <- read.csv(here("..", "data_clean", data_name), header = TRUE)

# load parameters
source(here("code", "params.R"))
dat <- dat.mock

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
  "Bstratum", "wt", "wt.2", "wt.subcohort", "race",
  "EthnicityHispanic","EthnicityNotreported", "EthnicityUnknown",
  "WhiteNonHispanic"
)] %>%
  replicate(length(assays), ., simplify = FALSE) %>%
  bind_rows()


dat.long.assay_value.names <- times
dat.long.assay_value <- as.data.frame(matrix(
  nrow = nrow(dat) * length(assays),
  ncol = length(dat.long.assay_value.names)
))
colnames(dat.long.assay_value) <- dat.long.assay_value.names

for (tt in seq_along(times)) {
  dat_mock_col_names <- paste(times[tt], assays, sep = "")
  dat.long.assay_value[, dat.long.assay_value.names[tt]] <- unlist(lapply(
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
      labels = c("Age < 65", "Age >= 65")
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
      levels = c(1, 0),
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
    ifelse(
      EthnicityHispanic == 1,
      "Hispanic or Latino",
      ifelse(
        EthnicityNotreported == 0 & EthnicityUnknown == 0,
        "Not Hispanic or Latino",
        "Not reported and unknown"
      )
    )
  ) %>% factor(
    levels = c("Hispanic or Latino", "Not Hispanic or Latino", "Not reported and unknown", "Others")
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

dat.long.twophase.sample$race <- as.factor(dat.long.twophase.sample$race)
dat.twophase.sample$race <- as.factor(dat.twophase.sample$race)

dat.long.twophase.sample$Ptid <- as.character(dat.long.twophase.sample$Ptid) 
dat.twophase.sample$Ptid <- as.character(dat.twophase.sample$Ptid) 

saveRDS(as.data.frame(dat.long.twophase.sample),
  file = here("data_clean", "long_twophase_data.rds")
)
saveRDS(as.data.frame(dat.twophase.sample),
  file = here("data_clean", "twophase_data.rds")
)

## "longer" dataset

