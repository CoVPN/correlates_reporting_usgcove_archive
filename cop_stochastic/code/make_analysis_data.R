#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

# load packages
library(here)
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")

# load data
data_name_amended <- paste0(str_remove(data_name, ".csv"),
                            "_with_riskscore.csv")
data_name_check <- file.exists(here("..", "data_clean", data_name_amended))
if (data_name_check) {
  full_data <- read.csv(here("..", "data_clean", data_name_amended))
} else {
  full_data <- read.csv(here("..", "data_clean", data_name))
}

# load parameters
source(here::here("code", "params.R"))
if (data_name_check) {
  covariates <- c(covariates, "risk_score")
}

# Generate the outcome and censoring indicator variables
for (time in times) {
  print(paste("Preparing estimation data for", time))

  # get event indicator, event time, and failure indicator for this timepoint
  this_event_ind <- event_ind[[time]]
  this_event_time <- event_time[[time]]
  this_tf_ind <- tf_ind[[time]]

  # subset markers for timepoint
  markers_at_time <- str_subset(markers, time)

  # make binary outcome variable for this timepoint
  outcome <- full_data %>%
    transmute(
      get(this_event_ind) == 1 & get(this_event_time) <= this_tf_ind
    ) %>%
    deframe() %>%
    as.numeric()

  # optionally define the loss to follow-up indicator
  if (adjust_censoring) {
    cens_ind <- full_data %>%
      transmute(
        1 - (get(this_event_ind) == 0 & get(this_event_time) < this_tf_ind)
      ) %>%
      deframe() %>%
      as.numeric()
  } else {
    cens_ind <- rep(1, length(outcome))
  }

  # make estimation-ready data
  data_for_est <- full_data %>%
    select(
      all_of(covariates),
      all_of(markers_at_time),
      twophase_ind[[time]],
      twophase_wts[[time]],
      twophase_group[[time]],
      "Trt",
      "Bserostatus"
    ) %>%
    rename(
      samp_ind = eval(twophase_ind[[time]]),
      samp_wts = eval(twophase_wts[[time]]),
      samp_grp = eval(twophase_group[[time]])
    ) %>%
    mutate(
      outcome = outcome,
      cens_ind = cens_ind,
    ) %>%
    # NOTE: used to have Perprotocol == 1 as well, still needed?
    filter(Bserostatus == 0, !is.na(samp_wts)) %>%
    select(-Bserostatus)

  # save estimation-ready data for this analysis
  saveRDS(
    object = data_for_est,
    file = here("data_clean", paste0("data_est_", time, ".rds"))
  )
}
