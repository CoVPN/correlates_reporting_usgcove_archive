#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(readr)
require(tidyverse)
library(here)

# Create fancy/tidy screen names for use in tables and figures
# @param avgs dataframe containing Screen, Learner, AUCs information as columns
# @return object containing tidy screen names
get_fancy_screen_names <- function(avgs) {
  return(avgs %>%
    mutate(fancyScreen = case_when(
      Screen == "screen_highcor_random_plus_exposure_All" ~ "highcor_random",
      Screen == "screen_glmnet_plus_exposure_All" ~ "glmnet",
      Screen == "screen_univariate_logistic_pval_plus_exposure_All" ~ "univar_logistic_pval",
      Screen == "screen_all_plus_exposure_All" ~ "all",
      TRUE ~ as.character(Screen)
    )))
}


# Drop seeds/fits that returned any error
# @param dat object containing all 10 fits (as lists) from the CV.Superlearner with folds and auc information
# @return object upon dropping any fit that returned an error
drop_seeds_with_error <- function(dat) {
  newdat <- vector(mode = "list", length = 1)
  j <- 1
  for (i in 1:10) {
    if (typeof(dat[[i]][1]) == "list") {
      newdat[[j]] <- dat[[i]]
      j <- j + 1
    }
  }
  newdat
}


# Convert SL object to SL results dataframe
# @param dat object containing all 10 fits (as lists) from the CV.Superlearner with folds and auc information
# @return dataframe containing CV-AUCs
convert_SLobject_to_Slresult_dataframe <- function(dat) {

  # Remove any iteration seeds that returned an error!
  newdat <- drop_seeds_with_error(dat)

  # if (is.null(newdat[[1]])) {
  #   return( read.csv("empty_df.csv", stringsAsFactors = FALSE) %>% select(-X) %>% as_tibble() )
  # }

  if (!is.null(newdat[[1]])) {
    as_tibble(do.call(rbind.data.frame, lapply(newdat, function(x) x$aucs))) %>%
      # filter(!is.na(ci_ll) | !is.na(ci_ul)) %>%    # drop learners that have NA for ci_ul or ci_ll for certain seeds!
      group_by(Learner, Screen) %>%
      summarise(AUC = mean(AUC), ci_ll = mean(ci_ll), ci_ul = mean(ci_ul), .groups = "drop") %>%
      arrange(-AUC) %>%
      mutate(
        AUCstr = paste0(format(round(AUC, 3), nsmall = 3), " [", format(round(ci_ll, 3), nsmall = 3), ", ", format(round(ci_ul, 3), nsmall = 3), "]"),
        Learner = as.character(Learner),
        Screen = as.character(Screen),
        LearnerScreen = paste(Learner, Screen)
      ) %>%
      get_fancy_screen_names() %>%
      rename(
        Screen_fromRun = Screen,
        Screen = fancyScreen
      )
  }
}



# Read in SL objects from folder, get AUCs in a dataframe
# @param data_file RDS file containing all 10 fits from the CV.Superlearner with folds and auc information
# @param trt string containing treatment arm (placebo or vaccine)
# @return dataframe containing CV-AUCs
readin_SLobjects_fromFolder <- function(data_file, trt) {
  readRDS(data_file) %>%
    convert_SLobject_to_Slresult_dataframe() %>%
    mutate(trt = trt)
}


# Read CV.SL object and save AUCs
data_file <- here("output", "cvsl_riskscore_cvaucs.rds")
risk_placebo_cvaucs <- readin_SLobjects_fromFolder(data_file, trt = "placebo")
save(risk_placebo_cvaucs, file = here("output", "cvsl_risk_placebo_cvaucs.rda"))
