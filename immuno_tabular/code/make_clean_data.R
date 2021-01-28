##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here(".."))       #
source(here::here("..", "_common.R"))            #
##################################################

# load packages and helper scripts
library(here)
library(tidyverse)
library(COVIDcorr)
source(here("code", "make_functions.R"))

# read in original data
data(dat.mock)
ds_o <- dat.mock

# The stratified random cohort for immunogenicity
ds_s <- ds_o %>%
  dplyr::filter(SubcohortInd == 1 & TwophasesampInd == 1 & Perprotocol == 1)

# Generate a long format dataset stacking the subgroups, so that
#    Bserostatus*Trt, Bserostatus*Trt*subgroup could be run together.
# ds_all: subgroup = "Total", by Bserostatus*Trt,
#         weight = inverse of sampling prob. of the 8 stratums.
#         subgroup = "MinorityInd", "HighRiskInd", "Agecat", "Sex",
#                    by Bserostatus*Trt*subgroup (respectively), weight = 1

# Email from Peter, 12/31/2020
# On page 36 below Table 6, it is said that "For purposes of analysis White is
# defined as Race=White and Ethnicity=Not Hispanic or Latino. All of the other
# Race subgroups are defined solely by the Race variable." Does it mean when
# plotting by race, I should omit the data of white Hispanic and white Latino?
# I think the answer is yes, White will always mean race==White &
# ethnicity==Not Hispanic or Latino This follows what Moderna did.

ds_long <- ds_s %>%
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
    AgeMinorC = paste(Age65C, MinorityC)
  ) %>%
  pivot_longer(
    cols = c(
      Age65C, HighRiskC, AgeRiskC, SexC, AgeSexC, ethnicityC,
      RaceEthC, MinorityC, AgeMinorC
    ),
    names_to = "subgroup", values_to = "subgroup_cat"
  ) %>%
  dplyr::filter(!subgroup_cat %in% c("", "Age < 65 ", "Age >= 65 ")) %>%
  dplyr::filter(!(subgroup == "ethnicityC" &
                  subgroup_cat == "Not reported and unknown")
               )

ds_all <- bind_rows(
  ds_long,
  ds_s %>% mutate(subgroup = "Total", subgroup_cat = "Total")
)

# SAP Section 9.1.1 Antibody Marker Data
# Truncation - LLOQ, ULOQ
ds1 <- bind_cols(
  ds_all %>% # Truncate at LLOQ and ULOQ
    mutate_at(bAb_v, setLOQ, lloq = bAb_lloq, uloq = bAb_uloq) %>%
    mutate_at(grep("50", c(pnAb_v, lnAb_v), value = T), setLOQ,
              lloq = nAb50_lloq) %>%
    mutate_at(grep("80", c(pnAb_v, lnAb_v), value = T), setLOQ,
              lloq = nAb80_lloq),

  ds_all %>% # Include the original data
    select_at(c(bAb_v, pnAb_v, lnAb_v)) %>%
    rename_all(function(x) paste0(x, "_org"))
)

# Responders
ds2 <- bind_cols(
  ds1,
  pmap(list(
    data = replicate(length(c(bAb, pnAb, lnAb)), ds1, simplify = FALSE),
    bl = paste0("B", c(bAb, pnAb, lnAb)),
    post = paste0("Day29", c(bAb, pnAb, lnAb)),
    lloq = list(bAb_lloq, bAb_lloq, nAb50_lloq, nAb80_lloq, nAb50_lloq)
  ),
  .f = setResponder, folds = c(2, 4), responderFR = 4
  ) %>%
  do.call(cbind, .),

  pmap(list(
    data = replicate(length(c(bAb, pnAb, lnAb)), ds1, simplify = F),
    bl = paste0("B", c(bAb, pnAb, lnAb)),
    post = paste0("Day57", c(bAb, pnAb, lnAb)),
    lloq = list(bAb_lloq, bAb_lloq, nAb50_lloq, nAb80_lloq, nAb50_lloq)
  ),
  .f = setResponder, folds = c(2, 4), responderFR = 4
  ) %>%
  do.call(cbind, .),
)

sub_lev <- c("Total", "Age65C", "HighRiskC", "AgeRiskC", "SexC", "AgeSexC",
             "ethnicityC", "RaceEthC", "MinorityC", "AgeMinorC")
sub_lb <- c("Total", "Age", "High Risk", "Age, Risk", "Sex", "Age, Sex",
            "Ethnicity", "Race", "Minority", "Age, Minority")
grp_lev <- ds2 %>%
  mutate(subgroup = factor(subgroup, levels = sub_lev)) %>%
  arrange(subgroup, subgroup_cat) %>%
  distinct(subgroup_cat) %>%
  pull(subgroup_cat)

ds <- ds2 %>%
  mutate(
    subgroup = factor(subgroup, levels = sub_lev, labels = sub_lb),
    Group = factor(subgroup_cat, levels = grp_lev),
    Baseline = ifelse(Bserostatus == 1, "Positive", "Negative")
  )

# Add labels
ep_lev <- c("Anti-Spike IgG", "Anti-RBD IgG", "Pseudo nAb ID50",
            "Pseudo nAb ID80", "Live Virus nAb ID50", "Live Virus nAb ID80")
resp_lb <- expand.grid(
  time = visits,
  endpoint = c(bAb, pnAb, lnAb),
  ind = c("FR2", "FR4", "Resp"),
  stringsAsFactors = F
) %>%
  mutate(
    Visit = case_when(time == "B" ~ "D1", time == "Day29" ~ "D29",
                      time == "Day57" ~ "D57"),
    Endpoint = case_when(
      endpoint == "bindSpike" ~ "Anti-Spike IgG",
      endpoint == "bindRBD" ~ "Anti-RBD IgG",
      endpoint == "pseudoneutid50" ~ "Pseudo nAb ID50",
      endpoint == "pseudoneutid80" ~ "Pseudo nAb ID80",
      endpoint == "liveneutmn50" ~ "Live Virus nAb MN50"
    ),
    Endpoint = factor(Endpoint, levels = ep_lev),
    ind.lb = case_when(ind == "FR2" ~ "2-Fold Rise",
                       ind == "FR4" ~ "4-Fold Rise",
                       ind == "Resp" ~ "Responder")
  )

resp_v <- grep("Resp|FR2|FR4", names(ds), value = TRUE)

ds_resp_l <- pivot_longer(ds, cols = all_of(resp_v),
                          names_to = "responder_cat",
                          values_to = "response") %>%
  inner_join(resp_lb %>% unite("responder_cat",
                               c(time, endpoint, ind), sep = ""),
             by = "responder_cat")

ds_mag_l <- pivot_longer(ds, cols = all_of(c(bAb_v, pnAb_v, lnAb_v)),
                         names_to = "responder_cat", values_to = "mag") %>%
  inner_join(resp_lb %>%
    distinct(time, endpoint, Visit, Endpoint) %>%
    unite("responder_cat", c(time, endpoint), sep = ""), by = "responder_cat")

save(ds, ds_resp_l, ds_mag_l, resp_v,
  file = here("data_clean", "ds_all.Rdata")
)
