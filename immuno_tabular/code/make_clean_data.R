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
dat.mock <- read.csv(here::here("../data_clean", data_name))

# The stratified random cohort for immunogenicity
ds_s <- dat.mock %>%
  dplyr::filter(SubcohortInd == 1 & TwophasesampInd == 1 & Perprotocol == 1) %>%
  dplyr::filter(!is.na(wt.subcohort)) %>%
  # select(-contains("liveneut")) %>%
  # The subgroup variables need to be character not factors
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
                        EventIndPrimaryD57 == 0 ~ "Per-protocol non-cases")
  ) 

# Generate a long format dataset stacking the subgroups

subgrp <- c(
  "Total" = "All participants", 
  "Age65C" = "Age",
  "HighRiskC" = "Risk for Severe Covid-19",
  "AgeRiskC" = "Age, Risk for Severe Covid-19",
  "SexC" = "Sex", 
  "AgeSexC" = "Age, sex",
  "ethnicityC" = "Hispanic or Latino ethnicity", 
  "RaceEthC" = "Race",
  "MinorityC" = "Underrepresented minority status",
  "AgeMinorC" = "Age, Underrepresented minority status"
)

ds_long <- ds_s %>%
  pivot_longer(
    cols = c(
      Age65C, HighRiskC, AgeRiskC, SexC, AgeSexC, ethnicityC, RaceEthC,
      MinorityC, AgeMinorC),
    names_to = "subgroup", values_to = "subgroup_cat") %>%
  mutate(subgroup = factor(subgrp[subgroup], levels = subgrp)) %>% 
  dplyr::filter(!is.na(subgroup_cat) & subgroup_cat != "White")

ds_all <- bind_rows(
  ds_long,
  ds_s %>% mutate(
    subgroup = factor("All participants", levels = subgrp),
    subgroup_cat = ""
  )
)

# SAP Section 9.1.1 Antibody Marker Data
# Step1: Truncation - LLOD, ULOQ
ds1 <- ds_all
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
    data = replicate(length(c(bAb, pnAb, lnAb))*post_n, ds1, simplify = FALSE),
    bl = as.vector(outer(rep("B", post_n), c(bAb, pnAb, lnAb), paste0)),
    post = as.vector(outer(post, c(bAb, pnAb, lnAb), paste0)),
    llod = llods[rep(c(bAb, pnAb, lnAb), each = post_n)]),
    .f = setResponder, folds = c(2, 4), responderFR = 4) %>%
    do.call(cbind, .),
  
  # % > 2lloq and 4lloq
  pmap(list(
    data = replicate(length(c(bAb_v, pnAb_v, lnAb_v)), ds1, simplify = FALSE),
    x = c(bAb_v, pnAb_v, lnAb_v),
    llod = llods[rep(c(bAb, pnAb, lnAb), each = (post_n + 1))]),
    .f = grtLLOD) %>%
    do.call(cbind, .)
)

# Step3: Delta
ds3 <- bind_cols(
  ds2,
  pmap(list(
    data = replicate(length(c(bAb, pnAb, lnAb)), ds2, simplify = FALSE),
    timepoints = replicate(length(c(bAb, pnAb, lnAb)), c("B", "Day29", "Day57"),
                           simplify = FALSE),
    marker = c(bAb, pnAb, lnAb)
  ),
  .f = setDelta
  ) %>%
    do.call(cbind, .)
)

grp_lev <- c("", "Age $<$ 65",  "Age $\\geq$ 65", "At-risk", "Not at-risk", 
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

ds <- mutate(ds3, Group = factor(subgroup_cat, levels = grp_lev))

resp_v <- setdiff(grep("Resp|FR2|FR4|2llod|4llod", names(ds), value = T),
                  grep("liveneut", names(ds), value = T)
                  )

ds_resp_l <- pivot_longer(ds,
                          cols = all_of(resp_v), 
                          names_to = "resp_cat",
                          values_to = "response") %>%
  inner_join(labels_all, by = "resp_cat")

mag_v <- setdiff(
  c(bAb_v, pnAb_v, lnAb_v, grep("DeltaDay", names(ds), value = T)), 
  grep("liveneut", names(ds), value = T)
  )

ds_mag_l <- pivot_longer(ds,
                         cols = all_of(mag_v),
                         names_to = "mag_cat", values_to = "mag") %>%
  inner_join(labels_all %>% 
               distinct(mag_cat, time, marker, Visit, Marker, label.short), 
             by = "mag_cat")

save(ds_long, ds_resp_l, ds_mag_l, labels_all,
     file = here::here("data_clean", "ds_all.Rdata")
)
