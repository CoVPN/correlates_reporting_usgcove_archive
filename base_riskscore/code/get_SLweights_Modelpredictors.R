#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(tidyverse)
library(here)
library(conflicted)
library(glmnet)
library(xgboost)
library(ranger)
conflict_prefer("filter", "dplyr")

load(file = here("output", "objects_for_running_SL.rda"))
rm(Y, X_riskVars, weights, maxVar)

# Get Superlearner weights
load(file = here("output", "sl_riskscore_slfits.rda"))
sl_weights <- sl_riskscore_slfits$coef %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Learner") %>%
  rename(`Weights` = ".") %>%
  arrange(-Weights)

sl_weights %>%
  mutate(Weight = format(round(Weights, 3), nsmall = 3)) %>%
  mutate(
    Screen = paste0("screen_", sapply(strsplit(Learner, "_screen_"), `[`, 2)),
    Learner = sapply(strsplit(Learner, "_screen"), `[`, 1)
  ) %>%
  select(Learner, Screen, Weight) %>%
  write.csv(here("output", "SL_weights.csv"))

top_models <- sl_weights %>%
  .$Learner

# Get predictors selected in the models with highest weights
for (i in seq_along(top_models)) {
  if(top_models[i] %in% c("SL.glm_screen_univariate_logistic_pval", 
                          "SL.glm.interaction_screen_highcor_random",
                          "SL.glm_screen_all",
                          "SL.glm_screen_glmnet",
                          "SL.glm_screen_highcor_random",
                          "SL.glm.interaction_screen_glmnet",
                          "SL.glm.interaction_screen_univariate_logistic_pval",
                          "SL.gam_screen_glmnet",
                          "SL.gam_screen_univariate_logistic_pval",
                          "SL.gam_screen_highcor_random")) {
    model <- sl_riskscore_slfits[["fitLibrary"]][[top_models[i]]]$object$coefficients %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Predictors") %>%
      rename(`Coefficient` = ".") %>%
      mutate(
        `Odds Ratio` = exp(`Coefficient`),
        Learner = top_models[i])
  }

  if (top_models[i] %in% c("SL.glmnet_screen_all")) {
    model <- coef(sl_riskscore_slfits[["fitLibrary"]][[top_models[i]]]$object) %>%
      as.matrix() %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Predictors") %>%
      rename(`Coefficient` = "1") %>%
      mutate(`Odds Ratio` = exp(`Coefficient`),
             Learner = top_models[i])
  }

  if (top_models[i] %in% c("SL.xgboost_screen_all")) {
    model <- xgboost::xgb.importance(model = sl_riskscore_slfits[["fitLibrary"]][[top_models[i]]]$object) %>%
      as.data.frame() %>%
      mutate(Learner = top_models[i])
  }

  if (top_models[i] %in% c("SL.ranger.imp_screen_all")) {
    model <- sl_riskscore_slfits[["fitLibrary"]][[top_models[i]]]$object$variable.importance %>%
      as.data.frame() %>%
      rename(Importance = ".") %>%
      tibble::rownames_to_column(var = "Predictors") %>%
      mutate(Learner = top_models[i])
  }

  if (top_models[i] == "SL.mean_screen_all")
	next

  if (i == 1) {
    all_models <- model
  } else {
    all_models <- bind_rows(all_models, model)
  }
}

options(scipen=999)

if(run_prod){
  all_models %>%
    left_join(sl_weights, by = "Learner") %>%
    mutate(
      Weight = format(round(Weights, 3), nsmall = 3),
      Coefficient = format(round(Coefficient, 3), nsmall = 3),
      `Odds Ratio` = format(round(`Odds Ratio`, 3), nsmall = 3),
      Importance = format(round(Importance, 3), nsmall = 3),
      Gain = format(round(Gain, 3), nsmall = 3),
      Cover = format(round(Cover, 3), nsmall = 3),
      Frequency = format(round(Frequency, 3), nsmall = 3),
    ) %>%
    mutate(
      Screen = paste0("screen_", sapply(strsplit(Learner, "_screen_"), `[`, 2)),
      Learner = sapply(strsplit(Learner, "_screen"), `[`, 1)
    ) %>%
    select(Learner, Screen, Weight, Predictors, Coefficient, `Odds Ratio`,
           Importance, Feature, Gain, Cover, Frequency) %>%
    write.csv(here("output", "SL_all_models_with_predictors.csv"))
}else{
  all_models %>%
    left_join(sl_weights, by = "Learner") %>%
    mutate(
      Weight = format(round(Weights, 3), nsmall = 3),
      Coefficient = format(round(Coefficient, 3), nsmall = 3),
      `Odds Ratio` = format(round(`Odds Ratio`, 3), nsmall = 3),
    ) %>%
    mutate(
      Screen = paste0("screen_", sapply(strsplit(Learner, "_screen_"), `[`, 2)),
      Learner = sapply(strsplit(Learner, "_screen"), `[`, 1)
    ) %>%
    select(Learner, Screen, Weight, Predictors, Coefficient, `Odds Ratio`) %>%
    write.csv(here("output", "SL_all_models_with_predictors.csv"))
}


rm(sl_riskscore_slfits)

