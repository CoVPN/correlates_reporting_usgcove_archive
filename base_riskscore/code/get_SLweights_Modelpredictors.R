#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(tidyverse)
library(here)
library(conflicted)
conflict_prefer("filter", "dplyr")

# Get Superlearner weights
load(file = here("output", "sl_riskscore_slfits.rda"))
sl_weights <- sl_riskscore_slfits$coef %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Learner") %>%
  rename(`Weights` = ".") %>%
  arrange(-Weights)

sl_weights %>%
  mutate(Weights = format(round(Weights, 3), nsmall = 3)) %>%
  mutate(
    Learner = str_remove(Learner, "_plus_exposure"),
    Learner = str_remove(Learner, "_All")
  ) %>%
  write.csv(here("output", "SL_weights.csv"))

top_models <- sl_weights %>%
  filter(Weights > 0.5) %>%
  .$Learner

# Get predictors selected in the models with highest weights
library(glmnet)

for (i in 1:length(top_models)) {
  model <- sl_riskscore_slfits[["fitLibrary"]][[top_models[i]]]$object$coefficients %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Predictors") %>%
    rename(`Coefficient` = ".") %>%
    mutate(
      `Odds Ratio` = exp(`Coefficient`),
      Learner = top_models[i]
    )

  if (i == 1) {
    all_models <- model
  } else {
    all_models <- bind_rows(all_models, model)
  }
}

all_models %>%
  left_join(sl_weights, by = "Learner") %>%
  select(Learner, Weights, Predictors, Coefficient, `Odds Ratio`) %>%
  mutate(
    `SL Weights` = format(round(Weights, 3), nsmall = 3),
    Coefficient = format(round(Coefficient, 3), nsmall = 3),
    `Odds Ratio` = format(round(`Odds Ratio`, 3), nsmall = 3)
  ) %>%
  mutate(
    Learner = str_remove(Learner, "_plus_exposure"),
    Learner = str_remove(Learner, "_All")
  ) %>%
  select(-Weights) %>%
  select(Learner, `SL Weights`, everything()) %>%
  write.csv(here("output", "SL_top_models_with_predictors.csv"))

rm(sl_riskscore_slfits)
