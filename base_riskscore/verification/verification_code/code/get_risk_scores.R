# Predict probability of outcome using SuperLearner	

#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here("..", "..", ".."))
source(here::here("..", "..", "..", "_common.R"))
#-----------------------------------------------

library(SuperLearner)
library(here)

sl_fit <- readRDS(here("output", "sl_fit.rds"))
newx <- readRDS(here("output", "newx.rds"))

all_pred <- predict(sl_fit, newdata = newx, onlySL = TRUE)
sl_pred <- all_pred$pred
risk_scores <- qlogis(sl_pred)

save(risk_scores, file = here("output", "risk_scores.rds"))