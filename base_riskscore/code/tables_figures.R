#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

## load libraries and source files #############################################
library(cvAUC)
library(conflicted)
library(tidyverse)
library(vimp)
library(kyotil)
library(grid)
library(gridExtra)
library(cowplot)
library(here)
conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
source(here("code", "utils.R"))
method <- "method.CC_nloglik" # since SuperLearner relies on this to be in GlobalEnv
ggplot2::theme_set(theme_cowplot())

load(file = here("output", "objects_for_running_SL.rda"))
rm(Y, X_riskVars, weights, maxVar)

load(file = here("output", "cvsl_risk_placebo_cvaucs.rda"))

######## Table of demographic variables used to derive the risk score ##########
dat <- read.csv(here::here("..", "data_clean", data_name)) %>%
  select(all_of(risk_vars))

dat %>%
  map(~ sum(is.na(.))) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Variable Name") %>%
  mutate(
    V1 = paste0(V1, "/", nrow(dat), " (", format(round((V1 / nrow(dat)) * 100, 1), nsmall = 1), "%)"),
    Definition = case_when(
      `Variable Name` == "Age" ~ "Age at enrollment in years, between 18 and 85",
      `Variable Name` == "Sex" ~ "Sex assigned at birth (1=female, 0=male)",
      `Variable Name` == "BMI" ~ "BMI at enrollment (kg/m^2)",
      `Variable Name` == "MinorityInd" ~ "Baseline covariate underrepresented minority status (1=minority, 0=non-minority)",
      `Variable Name` == "EthnicityHispanic" ~ "Indicator ethnicity = Hispanic (0 = Non-Hispanic)",
      `Variable Name` == "EthnicityNotreported" ~ "Indicator ethnicity = Not reported (0 = Non-Hispanic)",
      `Variable Name` == "EthnicityUnknown" ~ "Indicator ethnicity = Unknown (0 = Non-Hispanic)",
      `Variable Name` == "Black" ~ "Indicator race = Black (0 = White)",
      `Variable Name` == "Asian" ~ "Indicator race = Asian (0 = White)",
      `Variable Name` == "NatAmer" ~ "Indicator race = American Indian or Alaska Native (0 = White)",
      `Variable Name` == "PacIsl" ~ "Indicator race = Native Hawaiian or Other Pacific Islander (0 = White)",
      `Variable Name` == "WhiteNonHispanic" ~ "Indicator race = White or Caucasian (1 = White)",
      `Variable Name` == "Multiracial" ~ "Indicator race = Multiracial (0 = White)",
      `Variable Name` == "Other" ~ "Indicator race = Other (0 = White)",
      `Variable Name` == "Notreported" ~ "Indicator race = Not reported (0 = White)",
      `Variable Name` == "Unknown" ~ "Indicator race = unknown (0 = White)",
      `Variable Name` == "HighRiskInd" ~ "Baseline covariate high risk pre-existing condition (1=yes, 0=no)"
    ),
    Comments = ""
  ) %>%
  rename(`Total missing values` = V1) %>%
  select(`Variable Name`, Definition, `Total missing values`, Comments) %>%
  write.csv(here("output", "risk_vars.csv"))

######## learner-screens #######################################################
caption <- "All learner-screen combinations (28 in total) used as input to the superlearner."

if (run_prod) {
  tab <- risk_placebo_cvaucs %>%
    filter(!Learner %in% c("SL", "Discrete SL")) %>%
    select(Learner, Screen) %>%
    mutate(
      Screen = fct_relevel(Screen, c(
        "all", "glmnet", "univar_logistic_pval",
        "highcor_random"
      )),
      Learner = as.factor(Learner),
      Learner = fct_relevel(Learner, c(
        "SL.mean", "SL.glm", "SL.glm.interaction",
        "SL.glmnet", "SL.gam", 
        "SL.xgboost", "SL.ranger.imp"
      ))
    ) %>%
    arrange(Learner, Screen) %>%
    distinct(Learner, Screen) %>%
    rename("Screen*" = Screen)
} else {
  tab <- risk_placebo_cvaucs %>%
    filter(!Learner %in% c("SL", "Discrete SL")) %>%
    select(Learner, Screen) %>%
    mutate(
      Screen = fct_relevel(Screen, c(
        "all", "glmnet", "univar_logistic_pval",
        "highcor_random"
      )),
      Learner = as.factor(Learner),
      Learner = fct_relevel(Learner, c("SL.mean", "SL.glm"))
    ) %>%
    arrange(Learner, Screen) %>%
    distinct(Learner, Screen) %>%
    rename("Screen*" = Screen)
}

tab %>% write.csv(here("output", "learner-screens.csv"))

######## SLperformance-plac ####################################################
caption <- "Performance of Superlearner and all learner-screen combinations (CV-AUCs with 95\\% CIs) for risk score analyses using placebo group and EventIndPrimaryD57 as outcome. Constraint of np/20 is applied to all learners such that no more than 6 input variables were allowed in any model."

sl.perf <- risk_placebo_cvaucs %>%
  mutate(AUCstr = ifelse(AUC %in% tail(sort(AUC), 1), paste0(AUCstr, "*"), AUCstr)) %>%
  select(Learner, Screen, AUCstr)

sl.perf %>% write.csv(here("output", "SLperformance-plac.csv"))

################################################################################
# Forest plots for risk_placebo model, yd57 endpoint
options(bitmapType = "cairo")
if (run_prod) {
  png(file = here("figs", "risk_placebo_cvaucs.png"),
      width = 2000, height = 1100)
  top_learner <- make_forest_plot_prod(risk_placebo_cvaucs)
} else {
  png(file = here("figs", "risk_placebo_cvaucs.png"),
      width = 2000, height = 700)
  top_learner <- make_forest_plot_demo(risk_placebo_cvaucs)
}
grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot,
             ncol = 2)
dev.off()

################################################################################
# plot ROC curve and pred.Prob with SL, Discrete SL and top 2 best-performing individual Learners
top2_plac <- bind_rows(
  risk_placebo_cvaucs %>% arrange(-AUC) %>%
    filter(!Learner %in% c("SL", "Discrete SL")) %>%
    dplyr::slice(1:2),
  risk_placebo_cvaucs %>%
    filter(Learner == "SL"),
  risk_placebo_cvaucs %>%
    filter(Learner == "Discrete SL")
) %>%
  mutate(LearnerScreen = ifelse(Learner == "SL", "Super Learner",
                                ifelse(Learner == "Discrete SL", Learner,
                                       paste0(Learner, "_", Screen_fromRun))))

# Get cvsl fit and extract cv predictions
load(file = here("output", "cvsl_riskscore_cvfits.rda"))
pred <- get_cv_predictions(cv_fit = cvfits[[1]], cvaucDAT = top2_plac)

# plot ROC curve
options(bitmapType = "cairo")
png(file = here("figs", "ROCcurve_riskscore_plac.png"),
    width = 1000, height = 1000)
p1 <- plot_roc_curves(pred, cvaucDAT = top2_plac)
print(p1)
dev.off()

# plot pred prob plot
options(bitmapType = "cairo")
png(file = here("figs", "predProb_riskscore_plac.png"),
    width = 1100, height = 1400)
p2 <- plot_predicted_probabilities(pred)
print(p2)
dev.off()

# Use SuperLearner to generate risk scores!
load(file = here("output", "risk_placebo_ptids.rda"))
plac <- bind_cols(
  risk_placebo_ptids,
  pred %>% filter(Learner == "SL") %>% select(pred, AUCchar)
) %>%
  mutate(risk_score = log(pred / (1 - pred)))

write.csv(plac, here("output", "placebo_ptids_with_riskscores.csv"),
          row.names = FALSE)
rm(cvfits, pred, p1, p2)
