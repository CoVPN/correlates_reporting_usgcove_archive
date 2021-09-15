#-----------------------------------------------
# obligatory to append to the top of each script
here::i_am("base_riskscore/code/createRDAfiles_fromSLobjects.R")
renv::activate(project = here::here())

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("_common.R"))
#-----------------------------------------------

library(here)
library(readr)
require(tidyverse)

# Create fancy/tidy screen names for use in tables and figures
# @param avgs dataframe containing Screen, Learner, AUCs information as columns
# @return object containing tidy screen names
get_fancy_screen_names <- function(avgs) {
  return(avgs %>%
    mutate(fancyScreen = case_when(
      Screen == "screen_highcor_random" ~ "highcor_random",
      Screen == "screen_glmnet" ~ "glmnet",
      Screen == "screen_univariate_logistic_pval" ~ "univar_logistic_pval",
      Screen == "screen_all" ~ "all",
      Screen == "All" ~ "-",
      TRUE ~ as.character(Screen)
    )))
}


# Drop seeds/fits that returned any error
# @param dat object containing all 10 fits (as lists) from the CV.Superlearner with folds and auc information
# @return object upon dropping any fit that returned an error
drop_seeds_with_error <- function(dat) {
  newdat <- vector(mode = "list", length = 1)
  j <- 1
  for (i in seq_along(dat)) {
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
  # Remove any iteration seeds that returned an error (if any)!
  newdat <- drop_seeds_with_error(dat)
  if (!is.null(newdat[[1]])) {
    as_tibble(do.call(rbind.data.frame, lapply(newdat, function(x) x$aucs))) %>%
      group_by(Learner, Screen) %>%
      summarise(AUC = mean(AUC), se = sqrt(mean(se ^ 2)),
                .groups = "drop") %>%
      arrange(-AUC) %>%
      mutate(
        ci_ll = AUC - qnorm(0.975) * se, ci_ul = AUC + qnorm(0.975) * se,
        AUCstr = paste0(format(round(AUC, 3), nsmall = 3), " [",
                        format(round(ci_ll, 3), nsmall = 3), ", ",
                        format(round(ci_ul, 3), nsmall = 3), "]"),
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
data_file <- here("base_riskscore", "output", "cvsl_riskscore_cvaucs.rds")
risk_placebo_cvaucs <- readin_SLobjects_fromFolder(data_file, trt = "placebo")
save(risk_placebo_cvaucs, file = here("base_riskscore", "output", "cvsl_risk_placebo_cvaucs.rda"))
