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

print("Data preprocess")

# For immunogenicity characterization, complete ignore any information on cases
# vs. non-cases.  The goal is to characterize immunogenicity in the random
# subcohort, which is a stratified sample of enrolled participants. So,
# immunogenicity analysis is always done in ppts that meet all of the criteria.
dat.twophase.sample <- dat %>%
  filter(ph2.immuno == 1)
twophase_sample_id <- dat.twophase.sample$Ptid

if (study_name_code == "ENSEMBLE") {
  important.columns <- c("Ptid", "Trt", "MinorityInd", "HighRiskInd", "Age", "Sex",
    "Bserostatus", "Senior", "Bstratum", "wt.subcohort", 
    "race","EthnicityHispanic","EthnicityNotreported", 
    "EthnicityUnknown", "WhiteNonHispanic", "Country", "HIVinfection")
} else {
  important.columns <- c("Ptid", "Trt", "MinorityInd", "HighRiskInd", "Age", "Sex",
               "Bserostatus", "Senior", "Bstratum", "wt.subcohort", 
               "race","EthnicityHispanic","EthnicityNotreported", 
               "EthnicityUnknown", "WhiteNonHispanic")
}
## arrange the dataset in the long form, expand by assay types
## dat.long.subject_level is the subject level covariates;
## dat.long.assay_value is the long-form time variables that differ by the assay type
dat.long.subject_level <- dat[, important.columns] %>%
  replicate(length(assay_immuno), ., simplify = FALSE) %>%
  bind_rows()


dat.long.assay_value.names <- times
dat.long.assay_value <- as.data.frame(matrix(
  nrow = nrow(dat) * length(assay_immuno),
  ncol = length(dat.long.assay_value.names)
))
colnames(dat.long.assay_value) <- dat.long.assay_value.names

for (tt in seq_along(times)) {
  dat_mock_col_names <- paste(times[tt], assay_immuno, sep = "")
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

dat.long.assay_value$assay <- rep(assay_immuno, each = nrow(dat))

dat.long <- cbind(dat.long.subject_level, dat.long.assay_value)


## change the labels of the factors for plot labels
dat.long$Trt <- factor(dat.long$Trt, levels = c(0, 1), labels = trt.labels)
dat.long$Bserostatus <- factor(dat.long$Bserostatus,
  levels = c(0, 1),
  labels = bstatus.labels
)
dat.long$assay <- factor(dat.long$assay, levels = assay_immuno, labels = assay_immuno)

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
    factor(Senior,
      levels = c(0, 1),
      labels = paste0(c("Age < ", "Age >= "), age.cutoff)
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
    factor(paste0(Senior, HighRiskInd),
      levels = c("00", "01", "10", "11"),
      labels = c(
        paste0("Age < ", age.cutoff, " not at risk"),
        paste0("Age < ", age.cutoff, " at risk"),
        paste0("Age >= ", age.cutoff, " not at risk"),
        paste0("Age >= ", age.cutoff, " at risk")
      )
    )
  )

if (study_name_code == "COVE") {

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
      factor(paste0(Senior, Sex),
        levels = c("00", "01", "10", "11"),
        labels = c(
          paste0("Age < ", age.cutoff, " male"),
          paste0("Age < ", age.cutoff, " female"),
          paste0("Age >= ", age.cutoff, " male"),
          paste0("Age >= ", age.cutoff, " female")
        )
      )
    )

} else if (study_name_code == "ENSEMBLE") {
  
  dat.long.twophase.sample$sex_label <-
    with(
      dat.long.twophase.sample,
      factor(Sex,
             levels = c(0, 1, 2, 3),
             labels = c("Male", "Female", "Undifferentiated", "Unknown")
      )
    )
  
  dat.long.twophase.sample$age_sex_label <-
    with(
      dat.long.twophase.sample,
      factor(paste0(Senior, Sex),
             levels = c("00", "01", "02", "03", "10", "11", "12", "13"),
             labels = c(
               paste0("Age < ", age.cutoff, " male"),
               paste0("Age < ", age.cutoff, " female"),
               paste0("Age < ", age.cutoff, " undifferentiated"),
               paste0("Age < ", age.cutoff, " unknown"),
               paste0("Age >= ", age.cutoff, " male"),
               paste0("Age >= ", age.cutoff, " female"),
               paste0("Age >= ", age.cutoff, " undifferentiated"),
               paste0("Age >= ", age.cutoff, " unknown")
             )
      )
    )
  
  # Ignore undifferentiated participants
  dat.long.twophase.sample$sex_label <- ifelse(dat.long.twophase.sample$sex_label == 2,
                                               NA, dat.long.twophase.sample$sex_label)
  dat.long.twophase.sample$age_sex_label <- ifelse(endsWith(dat.long.twophase.sample$age_sex_label, 2),
                                                   NA, dat.long.twophase.sample$age_sex_label)
  
}

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
    factor(MinorityInd,
      levels = c(0, 1),
      labels = c("White Non-Hispanic", "Comm. of Color")
    )
  )

dat.long.twophase.sample$age_minority_label <-
  with(
    dat.long.twophase.sample,
    factor(paste0(Senior, MinorityInd),
      levels = c("01", "00", "11", "10"),
      labels = c(
        paste0("Age < ", age.cutoff, " Comm. of Color"),
        paste0("Age < ", age.cutoff, " White Non-Hispanic"),
        paste0("Age >= ", age.cutoff, " Comm. of Color"),
        paste0("Age >= ", age.cutoff, " White Non-Hispanic")
      )
    )
  )

if(study_name_code == "ENSEMBLE") {
  dat.long.twophase.sample$country_label <- factor(sapply(dat.long.twophase.sample$Country, function(x) {
    names(countries.ENSEMBLE)[countries.ENSEMBLE==x]
  }), levels = names(countries.ENSEMBLE))
  
  dat.long.twophase.sample$hiv_label <- factor(sapply(dat.long.twophase.sample$HIVinfection, function(x) {
    ifelse(x,
           "HIV Positive",
           "HIV Negative")
  }), levels=c("HIV Negative", "HIV Positive"))
}

dat.long.twophase.sample$race <- as.factor(dat.long.twophase.sample$race)
dat.twophase.sample$race <- as.factor(dat.twophase.sample$race)

dat.long.twophase.sample$Ptid <- as.character(dat.long.twophase.sample$Ptid) 
dat.twophase.sample$Ptid <- as.character(dat.twophase.sample$Ptid) 


dat.long.twophase.sample <- filter(dat.long.twophase.sample, assay %in% assay_immuno)


saveRDS(as.data.frame(dat.long.twophase.sample),
  file = here("data_clean", "long_twophase_data.rds")
)
saveRDS(as.data.frame(dat.twophase.sample),
  file = here("data_clean", "twophase_data.rds")
)
