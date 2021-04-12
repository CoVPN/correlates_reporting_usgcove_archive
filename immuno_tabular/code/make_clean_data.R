##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here("..")) #
source(here::here("..", "_common.R")) #
##################################################

# Immunogenicity Tables

# Reload clean_data
base::load(here::here("data_clean", "params.Rdata"))

library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# Read in original data
dat.mock <- read.csv(here::here("..", "data_clean", data_name))

# The stratified random cohort for immunogenicity
ds_s <- dat.mock %>%
  dplyr::filter(Perprotocol == 1 & EventTimePrimaryD57 >= 7) %>%
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

    Case = case_when(Perprotocol == 1 & EventIndPrimaryD29 == 1 & 
                       EventIndPrimaryD57 == 1 ~ "Cases",
                     TRUE ~ "Non-Cases"),
    Case2 = case_when(EventIndPrimaryD29 == 1 & EventIndPrimaryD57 == 0 ~ 
                        "Per-protocol Intercurrent cases",
                      Perprotocol == 1 & EventIndPrimaryD29 == 1 & 
                        EventIndPrimaryD57 == 1 ~ "Per-protocol cases",
                      Perprotocol == 1 & EventIndPrimaryD29 == 0 & 
                        EventIndPrimaryD57 == 0 ~ "Per-protocol non-cases"),
    AgeRisk1 = ifelse(Age65C=="Age $<$ 65", AgeRiskC, NA),
    AgeRisk2 = ifelse(Age65C=="Age $\\geq$ 65", AgeRiskC, NA),
    All = "All participants",
    randomset = (SubcohortInd == 1 & TwophasesampInd == 1 &
                   Perprotocol == 1 & !is.na(wt.subcohort))
  ) 


# SAP Section 9.1.1 Antibody Marker Data
# Step1: Truncation - LLOD, ULOQ
ds1 <- ds_s
for(x in names(labels.assays.short)) {
  ds1 <- mutate_at(ds1, grep(x, names(ds1), value=T), 
                   setLOD, llod = llods[x], uloq = uloqs[x])
} 

# Step2: Responders
# Post baseline visits
post <- names(labels.time)[!grepl("B|Delta", names(labels.time), fixed = F)]
post_n <- length(post)

ds2 <- bind_cols(
  # Original ds
  ds1,
  # Responses post baseline
  pmap(list(
    data = replicate(length(assays)*post_n, ds1, simplify = FALSE),
    bl = as.vector(outer(rep("B", post_n), assays, paste0)),
    post = as.vector(outer(post, assays, paste0)),
    llod = llods[rep(assays, each = post_n)]),
    .f = setResponder, folds = c(2, 4), responderFR = 4) %>%
    do.call(cbind, .),
  
  # % > 2lloq and 4lloq
  pmap(list(
    data = replicate(length(assays_col), ds1, simplify = FALSE),
    x = assays_col,
    llod = llods[rep(assays, each = (post_n + 1))]),
    .f = grtLLOD) %>%
    do.call(cbind, .)
)

# Step3: Delta
ds <- bind_cols(
  ds2 %>% select(!contains("Delta")),
  pmap(list(
    data = replicate(length(assays), ds2, simplify = FALSE),
    timepoints = replicate(length(assays), c("B", "Day29", "Day57"),
                           simplify = FALSE),
    marker = assays
  ),
  .f = setDelta
  ) %>%
    do.call(cbind, .)
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
  MinorityC = "Underrepresented minority status",
  AgeMinorC = "Age, Underrepresented minority status"
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
