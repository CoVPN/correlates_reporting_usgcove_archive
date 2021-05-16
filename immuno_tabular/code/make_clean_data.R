##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here("..")) #
source(here::here("..", "_common.R")) #
##################################################

# Immunogenicity Tables

# Reload clean_data
base::load(here::here("data_clean", "params.Rdata"))
source(here::here("code", "make_functions.R"))

library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# Read in original data
dat <- read.csv(here::here("..", "data_clean", data_name))

# The stratified random cohort for immunogenicity
ds_s <- dat %>%
  dplyr::filter(!is.na(wt.subcohort)) %>%
  mutate(
    raceC = as.character(race),
    ethnicityC = case_when(EthnicityHispanic==1 ~ "Hispanic or Latino",
                           EthnicityHispanic==0 & EthnicityNotreported==0 & 
                             EthnicityUnknown==0 ~ "Not Hispanic or Latino",
                           EthnicityNotreported==1 | 
                             EthnicityUnknown==1 ~ "Not reported and unknown "),
    RaceEthC = case_when(
      WhiteNonHispanic==1 ~ "White Non-Hispanic ",
      TRUE ~ raceC
    ),
    MinorityC = case_when(
      WhiteNonHispanic == 0 ~ "Communities of Color",
      WhiteNonHispanic == 1 ~ "White Non-Hispanic"
    ),
    HighRiskC = ifelse(HighRiskInd == 1, "At-risk", "Not at-risk"),
    Age65C = ifelse(age.geq.65 == 1, "Age $\\geq$ 65", "Age $<$ 65"),
    SexC = ifelse(Sex == 1, "Female", "Male"),
    AgeRiskC = paste(Age65C, HighRiskC),
    AgeSexC = paste(Age65C, SexC),
    AgeMinorC = ifelse(is.na(MinorityC), NA, paste(Age65C, MinorityC)),
    `Baseline SARS-CoV-2` = factor(ifelse(Bserostatus == 1, "Positive", "Negative"),
                      levels = c("Negative", "Positive")
    ),
    Arm = factor(ifelse(Trt == 1, "Vaccine", "Placebo"), 
                levels = c("Vaccine", "Placebo")),
    AgeRisk1 = ifelse(Age65C=="Age $<$ 65", AgeRiskC, NA),
    AgeRisk2 = ifelse(Age65C=="Age $\\geq$ 65", AgeRiskC, NA),
    All = "All participants",
    randomset = (SubcohortInd == 1 & TwophasesampIndD57 == 1 & EarlyendpointD57==0)
  ) 

# Step2: Responders
# Post baseline visits
post <- names(labels.time)[!grepl("B|Delta", names(labels.time), fixed = F)]
post_n <- length(post)

ds <- bind_cols(
  # Original ds
  ds_s,
  # Responses post baseline
  pmap(list(
    data = replicate(length(assays)*post_n, ds_s, simplify = FALSE),
    bl = as.vector(outer(rep("B", post_n), assays, paste0)),
    post = as.vector(outer(post, assays, paste0)),
    cutoff = lloqs[rep(assays, each = post_n)]),
    .f = setResponder, folds = c(2, 4), responderFR = 4) %>%
    bind_cols(),
  
  # % > 2lloq and 4lloq
  pmap(list(
    data = replicate(length(assays_col), ds_s, simplify = FALSE),
    x = assays_col,
    cutoff.name="lloq",
    cutoff = lloqs[rep(assays, each = (post_n + 1))]),
    .f = grtLL) %>%
    bind_cols()
)

subgrp <- c(
  All = "All participants", 
  Age65C = "Age",
  BMI="BMI",
  HighRiskC = "Risk for Severe Covid-19",
  AgeRiskC = "Age, Risk for Severe Covid-19",
  AgeRisk1 = "Age < 65, Risk for Severe Covid-19",
  AgeRisk2 = "Age >= 65, Risk for Severe Covid-19",
  SexC = "Sex", 
  AgeSexC = "Age, sex",
  ethnicityC = "Hispanic or Latino ethnicity", 
  RaceEthC = "Race",
  MinorityC = "Communities of color",
  AgeMinorC = "Age, Communities of color"
)

grplev <- c("", "Age $<$ 65",  "Age $\\geq$ 65", "At-risk", "Not at-risk", 
            "Age $<$ 65 At-risk", "Age $<$ 65 Not at-risk", 
            "Age $\\geq$ 65 At-risk", "Age $\\geq$ 65 Not at-risk", "Male", "Female", 
            "Age $<$ 65 Female", "Age $<$ 65 Male", 
            "Age $\\geq$ 65 Male", "Age $\\geq$ 65 Female", 
            "Hispanic or Latino", "Not Hispanic or Latino", "Not reported and unknown ", 
            "White Non-Hispanic ", "Black or African American", "Asian", 
            "American Indian or Alaska Native", 
            "Native Hawaiian or Other Pacific Islander", 
            "Multiracial", "Other", "Not reported and unknown",  
            "Communities of Color", "White Non-Hispanic",
            "Age $<$ 65 Communities of Color", "Age $<$ 65 White Non-Hispanic",  
            "Age $\\geq$ 65 Communities of Color", "Age $\\geq$ 65 White Non-Hispanic")
names(grplev) <- c("All participants", grplev[-1])

save(ds, assays, assays_col, labels_all, subgrp, grplev, tlf, 
     file = here::here("data_clean", "ds_all.Rdata"))

