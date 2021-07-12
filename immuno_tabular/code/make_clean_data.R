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

if (subset_variable != "None") {
  dat <- dat %>% dplyr::filter(!!as.name(subset_variable)==subset_value)
  }

# The stratified random cohort for immunogenicity

ds_s <- dat %>%
  dplyr::filter(ph1.immuno) %>%
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
    AgeC = ifelse(Senior == 1, labels.age[2], labels.age[1]),
    SexC = ifelse(Sex == 1, "Female", "Male"),
    AgeRiskC = paste(AgeC, HighRiskC),
    AgeSexC = paste(AgeC, SexC),
    AgeMinorC = ifelse(is.na(MinorityC), NA, paste(AgeC, MinorityC)),
    `Baseline SARS-CoV-2` = factor(ifelse(Bserostatus == 1, "Positive", "Negative"),
                      levels = c("Negative", "Positive")
    ),
    Arm = factor(ifelse(Trt == 1, "Vaccine", "Placebo"), 
                levels = c("Vaccine", "Placebo")),
    AgeRisk1 = ifelse(AgeC==labels.age[1], AgeRiskC, NA),
    AgeRisk2 = ifelse(AgeC==labels.age[2], AgeRiskC, NA),
    All = "All participants"
  ) 

if(study_name_code=="ENSEMBLE"){
  ds_s <- ds_s %>% 
    mutate(CountryC=labels.countries.ENSEMBLE[Country+1],
           RegionC=labels.regions.ENSEMBLE[Region+1],
           URMC=case_when(URMforsubcohortsampling == 1 ~ "URM",
                          URMforsubcohortsampling == 0 ~ "Non-URM"),
           AgeURM=paste(AgeC, URMC))
}

# Step2: Responders, % >=2FR, % >=4FR, % >=2lloq, % >=4lloq
# Post baseline visits
ds <- getResponder(ds_s, cutoff.name=cutoff.name, times=grep("Day", times, value=T), 
                   assays=assays, pos.cutoffs=pos.cutoffs)

subgrp <- c(
  All = "All participants", 
  AgeC = "Age",
  BMI="BMI",
  HighRiskC = "Risk for Severe Covid-19",
  AgeRiskC = "Age, Risk for Severe Covid-19",
  AgeRisk1 = paste0(labels.age[1], ", Risk for Severe Covid-19"),
  AgeRisk2 = paste0(labels.age[2], ", Risk for Severe Covid-19"),
  SexC = "Sex", 
  AgeSexC = "Age, sex",
  ethnicityC = "Hispanic or Latino ethnicity", 
  RaceEthC = "Race",
  MinorityC = "Communities of color",
  AgeMinorC = "Age, Communities of color",
  URMC = "Underrepresented Minority Status in the U.S.",
  AgeURM = "Age, Underrepresented Minority Status in the U.S.",
  CountryC = "Country"
)

grplev <- c("", labels.age, "At-risk", "Not at-risk", 
            paste(labels.age[1], c("At-risk", "Not at-risk")),
            paste(labels.age[2], c("At-risk", "Not at-risk")),
            "Male", "Female", 
            paste(labels.age[1], c("Female", "Male")),
            paste(labels.age[2], c("Female", "Male")),
            "Hispanic or Latino", "Not Hispanic or Latino", "Not reported and unknown ", 
            "White Non-Hispanic ", "Black or African American", "Asian", 
            "American Indian or Alaska Native", 
            "Native Hawaiian or Other Pacific Islander", 
            "Multiracial", "Other", "Not reported and unknown",  
            labels.minor, 
            paste(labels.age[1], labels.minor),  
            paste(labels.age[2], labels.minor),
            labels.countries.ENSEMBLE)

names(grplev) <- c("All participants", grplev[-1])

save(ds, assays, assays_col, labels_all, subgrp, grplev, tlf, 
     file = here::here("data_clean", "ds_all.Rdata"))

