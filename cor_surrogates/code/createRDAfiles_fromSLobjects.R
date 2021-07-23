#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
#-----------------------------------------------

## ----createRDAfiles_fromSLobjects---------
library(readr)
require(tidyverse)
library(here)
source(here("code", "utils.R"))

# Create fancy/tidy screen names for use in tables and figures
# @param avgs dataframe containing Screen, Learner, AUCs information as columns
# @return object containing tidy screen names
get_fancy_screen_names <- function(avgs){
  return(avgs %>%
           mutate(fancyScreen = case_when(Screen == "screen_highcor_random" ~ "highcor_random",
                                          Screen == "screen_glmnet" ~ "glmnet",
                                          Screen == "screen_univariate_logistic_pval" ~ "univar_logistic_pval",
                                          Screen == "screen_all" ~ "all",
                                          TRUE ~ as.character(Screen)))
         )
}

# Drop seeds/fits that returned any error
# @param dat object containing all 10 fits (as lists) from the CV.Superlearner with folds and auc information
# @return object upon dropping any fit that returned an error
drop_seeds_with_error <- function(dat){
  newdat <- vector(mode = "list", length = 1)
  j = 1
  for (i in 1:length(dat)){
    if( typeof(dat[[i]][1]) == "list" ){
      newdat[[j]] = dat[[i]]
      j = j + 1
    }
  }
  newdat
}


# Convert SL object to SL results dataframe
# @param dat object containing all 10 fits (as lists) from the CV.Superlearner with folds and auc information
# @return dataframe containing CV-AUCs
convert_SLobject_to_Slresult_dataframe <- function(dat) {

  # Remove any iteration seeds that returned an error!
  newdat = drop_seeds_with_error(dat)
    
  if (is.null(newdat[[1]])) {
    return( read.csv("empty_df.csv", stringsAsFactors = FALSE) %>% select(-X) %>% as_tibble() )
  }

  if (!is.null(newdat[[1]])) {
    as_tibble(do.call(rbind.data.frame, lapply(newdat, function(x) x))) %>%
      filter(!is.na(ci_ll) | !is.na(ci_ul)) %>%    # drop learners that have NA for ci_ul or ci_ll for certain seeds!
      group_by(Learner, Screen) %>%
      summarize(AUC = mean(AUC), ci_ll = mean(ci_ll), ci_ul = mean(ci_ul), .groups = 'drop') %>%
      ungroup()  %>%
      arrange(-AUC) %>%
      mutate(AUCstr = paste0(format(round(AUC, 3), nsmall=3), " [", format(round(ci_ll, 3), nsmall=3), ", ", format(round(ci_ul, 3), nsmall=3), "]"),
             Learner = as.character(Learner),
             Screen = as.character(Screen),
             LearnerScreen = paste(Learner, Screen)) %>%
      get_fancy_screen_names() %>%
      rename(Screen_fromRun = Screen,
             Screen = fancyScreen)
  }
}



# Read in SL objects from folder, get AUCs in a dataframe
# @param data_file RDS file containing all 10 fits from the CV.Superlearner with folds and auc information
# @param trt string containing treatment arm (placebo or vaccine)
# @return dataframe containing CV-AUCs
readin_SLobjects_fromFolder <- function(data_path, file_pattern, endpoint, trt){
  dir(data_path, pattern = file_pattern) %>%
    tibble(file = .) %>%
    mutate(listdat = lapply(paste0(data_path, "/", file), readRDS)) %>% 
    mutate(data = map(listdat, convert_SLobject_to_Slresult_dataframe)) %>%
    select(file, data) %>%
    unnest(data) %>%
    mutate(endpoint = endpoint,
           trt = trt)
}


# Read CV.SL object and save relevant columns as dataframe
# For vaccine, yd57 endpoint
data_folder <- here("output")
cvaucs_vacc <- readin_SLobjects_fromFolder(data_folder, file_pattern = "*.rds", endpoint = "EventIndPrimaryD57", trt = "vaccine") %>%
  mutate(varset = str_replace(file, "CVSLaucs_vacc_EventIndPrimaryD57_", ""),
         varset = str_replace(varset, "_varset", ""),
         varset = str_replace(varset, ".rds", ""),
         varsetNo = as.numeric(sapply(strsplit(varset, "_"), `[`, 1))) %>%
  arrange(varsetNo)

save(cvaucs_vacc, file = here("output", "cvaucs_vacc_EventIndPrimaryD57.rda"))

