## Written by Bhavesh Borate
## For the RSV project
## dated 11 Dec 2020

## ----createRDAfiles_fromSLobjects---------
library(readr)
require(tidyverse)
library(here)
source(paste0(here(), "/correlates_report/SL_estimated_optimal_surrogate/cluster_code/utils.R"))

get_fancy_screen_names <- function(avgs){
  return(avgs %>%
           mutate(fancyScreen = case_when(Screen == "screen_highcor_random" ~ "highcor_random",
                                          Screen == "screen_glmnet" ~ "glmnet",
                                          Screen == "screen_univariate_logistic_pval" ~ "univar_logistic_pval",
                                          Screen == "screen_all" ~ "all",
                                          TRUE ~ as.character(Screen)))
         )
}



get_all_aucs_lst <- function(sl_fit_lst) {
  weights = rep(1, length(sl_fit_lst$fit$Y))
  # get the CV-R^2 of the SuperLearner predictions
  if (is.null(sl_fit_lst)) {
    return(NA)
  } else {
    sl_auc <- cv_auc(preds = sl_fit_lst$fit$SL.predict, Y = sl_fit_lst$fit$Y, folds = sl_fit_lst$fit$folds, weights = weights)
    out <- data.frame(Learner="SL", Screen="All", AUC = sl_auc$auc, ci_ll = sl_auc$ci[1], ci_ul=sl_auc$ci[2])
    
    # Get the CV-auc of the Discrete SuperLearner predictions
    discrete_sl_auc <- cv_auc(preds = sl_fit_lst$fit$discreteSL.predict, Y = sl_fit_lst$fit$Y, folds = sl_fit_lst$fit$folds, weights = weights)
    out <- rbind(out, data.frame(Learner="Discrete SL", Screen="All", AUC = discrete_sl_auc$auc, ci_ll = discrete_sl_auc$ci[1], ci_ul = discrete_sl_auc$ci[2]))
    
    # Get the cvauc of the individual learners in the library
    get_individual_auc <- function(sl_fit, col, weights) {
      if(any(is.na(sl_fit$library.predict[, col]))) return(NULL)
      alg_auc <- cv_auc(preds = sl_fit$library.predict[, col], Y = sl_fit$Y, folds = sl_fit$folds, weights = weights)
      ## get the regexp object
      alg_screen_string <- strsplit(colnames(sl_fit$library.predict)[col], "_", fixed = TRUE)[[1]]
      alg <- tail(alg_screen_string[grepl(".", alg_screen_string, fixed = TRUE)], n = 1)
      screen <- paste0(alg_screen_string[!grepl(alg, alg_screen_string, fixed = TRUE)], collapse = "_")
      data.frame(Learner = alg, Screen = screen, AUC = alg_auc$auc, ci_ll = alg_auc$ci[1], ci_ul = alg_auc$ci[2])
    }
    other_aucs <- plyr::ldply(1:ncol(sl_fit_lst$fit$library.predict), function(x) get_individual_auc(sl_fit_lst$fit, x, weights))
    return(rbind(out, other_aucs))
  }
}


drop_seeds_with_error <- function(dat){
  newdat <- vector(mode = "list", length = 1)
  j = 1
  for (i in 1:10){
    if( typeof(dat[[i]][1]) == "list" ){
      newdat[[j]] = dat[[i]]
      j = j + 1
    }
  }
  newdat
}

convert_SLobject_to_Slresult_dataframe <- function(dat) {

  # Remove any iteration seeds that returned an error!
  newdat = drop_seeds_with_error(dat)
    
  if (is.null(newdat[[1]])) {
    return( read.csv("empty_df.csv", stringsAsFactors = FALSE) %>% select(-X) %>% as_tibble() )
  }

  if (!is.null(newdat[[1]])) {
    as_tibble(do.call(rbind.data.frame, lapply(newdat, function(x) x$cvaucs$aucs))) %>%
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
cvaucs_d57_vacc <- readin_SLobjects_fromFolder(data_folder, file_pattern = "*.rds", endpoint = "EventIndPrimaryD57", trt = "vaccine") %>%
  mutate(file = str_replace(file, "CVSLfits_vacc_cvaucs_EventIndPrimaryD57_", ""),
         #file = str_replace(file, "slfits_", ""),
         #file = str_replace(file, "_y2_placebo", ""),
         file = str_replace(file, ".rds", "")) 
# %>%
#   separate(file, c("slrun", "varset"), "_vacc_") 

save(cvaucs_d57_vacc, file = here("output", "cvaucs_d57_vacc.rda"))

