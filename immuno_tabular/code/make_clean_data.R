##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here("..")) #
source(here::here("..", "_common.R")) #
##################################################

library(tidyverse)
# DB: Scheduled for deletion
# library(COVIDcorr)
load(here::here("data_clean", "params.Rdata"))

# DB: Scheduled for deletion
# # Read in original data
# data(dat.mock)
dat.mock <- read.csv(here("..", "data_raw", data_name))

ds_o <- dat.mock
# The stratified random cohort for immunogenicity
ds_s <- ds_o %>%
  dplyr::filter(SubcohortInd == 1 & TwophasesampInd == 1 & Perprotocol == 1) %>%
  # dplyr::filter(SubcohortInd == 1 & Perprotocol == 1) %>%
  mutate(raceC = as.character(race), ethnicityC = as.character(ethnicity)) %>%
  mutate(
    RaceEthC = case_when(
      raceC == "White" & ethnicityC == "Not Hispanic or Latino" ~
      "White Non-Hispanic",
      raceC == "White" & ethnicityC == "Hispanic or Latino" ~ "",
      TRUE ~ raceC
    ),
    MinorityC = case_when(
      MinorityInd == 1 ~ "Communities of Color",
      MinorityInd == 0 & raceC == "White" &
        ethnicityC == "Not Hispanic or Latino" ~ "White Non-Hispanic",
      TRUE ~ ""
    ),
    HighRiskC = ifelse(HighRiskInd == 1, "At-risk", "Not at-risk"),
    Age65C = ifelse(age.geq.65 == 1, "Age >= 65", "Age < 65"),
    SexC = ifelse(Sex == 1, "Female", "Male"),
    AgeRiskC = paste(Age65C, HighRiskC),
    AgeSexC = paste(Age65C, SexC),
    AgeMinorC = paste(Age65C, MinorityC),
    Baseline = factor(ifelse(Bserostatus == 1, "Positive", "Negative"),
      levels = c("Negative", "Positive")
    ),
    Rx = factor(ifelse(Trt == 1, "Vaccine", "Placebo"),
      levels = c("Vaccine", "Placebo")
    )
  )

# Generate a long format dataset stacking the subgroups, so that
# Bserostatus*Trt, Bserostatus*Trt*subgroup, etc., could be run together.
# ds_all: subgroup = "Total", by Bserostatus*Trt, weight=inverse of sampling
#                    probability of the 8 stratums.
#         subgroup = "MinorityInd", "HighRiskInd", "Agecat", "Sex", by
#                    Bserostatus*Trt*subgroup(respectively), weight=1

subgrp <- c(
  "Total" = "All participants", "Age65C" = "Age",
  "HighRiskC" = "Risk for Severe Covid-19",
  "AgeRiskC" = "Age, Risk for Severe Covid-19",
  "SexC" = "Sex", "AgeSexC" = "Age, sex",
  "ethnicityC" = "Hispanic or Latino ethnicity", "RaceEthC" = "Race",
  "MinorityC" = "Race and ethnic group",
  "AgeMinorC" = "Age, Race and ethnic group"
)

ds_long <- ds_s %>%
  pivot_longer(
    cols = c(
      Age65C, HighRiskC, AgeRiskC, SexC, AgeSexC, ethnicityC, RaceEthC,
      MinorityC, AgeMinorC
    ),
    names_to = "subgroup", values_to = "subgroup_cat"
  ) %>%
  mutate(subgroup = factor(subgrp[subgroup], levels = subgrp)) %>%
  dplyr::filter(!subgroup_cat %in% c("", "Age < 65 ", "Age >= 65 ")) %>%
  dplyr::filter(!(subgroup == "ethnicityC" &
    subgroup_cat == "Not reported and unknown"))

ds_all <- bind_rows(
  ds_long,
  ds_s %>% mutate(
    subgroup = factor("All participants", levels = subgrp),
    subgroup_cat = "All participants"
  )
)

# SAP Section 9.1.1 Antibody Marker Data
# Truncation - LLOQ, ULOQ
ds1 <- bind_cols(
  ds_all %>% # Truncate at LLOQ and ULOQ
    mutate_at(bAb_v, setLOQ, lloq = bAb_lloq, uloq = bAb_uloq) %>%
    mutate_at(grep("50", c(pnAb_v, lnAb_v), value = T), setLOQ,
      lloq = .getLLOQ("50")
    ) %>%
    mutate_at(grep("80", c(pnAb_v, lnAb_v), value = T), setLOQ,
      lloq = .getLLOQ("80")
    ),

  ds_all %>% # Include the original data
    select_at(c(bAb_v, pnAb_v, lnAb_v)) %>%
    rename_all(function(x) paste0(x, "_org"))
)

# Responders
ds2 <- bind_cols(
  # Original ds
  ds1,
  # Responses at Day 29
  pmap(list(
    data = replicate(length(c(bAb, pnAb, lnAb)), ds1, simplify = FALSE),
    bl = paste0("B", c(bAb, pnAb, lnAb)),
    post = paste0("Day29", c(bAb, pnAb, lnAb)),
    lloq = lapply(c(bAb, pnAb, lnAb), .getLLOQ)
  ),
  .f = setResponder, folds = c(2, 4), responderFR = 4
  ) %>%
    do.call(cbind, .),

  # Responses at Day 57
  pmap(list(
    data = replicate(length(c(bAb, pnAb, lnAb)), ds1, simplify = FALSE),
    bl = paste0("B", c(bAb, pnAb, lnAb)),
    post = paste0("Day57", c(bAb, pnAb, lnAb)),
    lloq = lapply(c(bAb, pnAb, lnAb), .getLLOQ)
  ),
  .f = setResponder, folds = c(2, 4), responderFR = 4
  ) %>%
    do.call(cbind, .)
)

# Delta
ds3 <- bind_cols(
  ds2,
  pmap(list(
    data = replicate(length(c(bAb, pnAb, lnAb)), ds2, simplify = FALSE),
    t1 = "B",
    t2 = "Day57",
    endpoint = c(bAb, pnAb, lnAb)
  ),
  .f = setDelta
  ) %>%
    do.call(cbind, .)
)

grp_lev <- ds2 %>%
  arrange(subgroup, subgroup_cat) %>%
  distinct(subgroup_cat) %>%
  pull(subgroup_cat)

ds <- ds2 %>%
  mutate(Group = factor(subgroup_cat, levels = grp_lev))

resp_v <- grep("Resp|FR2|FR4", names(ds), value = T)

ds_resp_l <- pivot_longer(ds,
  cols = all_of(resp_v), names_to = "resp_cat",
  values_to = "response"
) %>%
  inner_join(labels_all, by = "resp_cat")

ds_mag_l <- pivot_longer(ds,
  cols = all_of(c(bAb_v, pnAb_v, lnAb_v)),
  names_to = "mag_cat", values_to = "mag"
) %>%
  inner_join(labels_all %>% distinct(
    mag_cat, time, endpoint, Visit, Endpoint,
    label.long, label.short
  ), by = "mag_cat")


save(ds, ds_s, ds_long, ds_resp_l, ds_mag_l, labels_all,
  file = here::here("data_clean", "ds_all.Rdata")
)
