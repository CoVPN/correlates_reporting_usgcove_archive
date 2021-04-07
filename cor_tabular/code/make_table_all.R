##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here("..")) #
source(here::here("..", "_common.R")) #
##################################################

library(survey)
library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
# For stratum with 1 ppt
options(survey.lonely.psu="adjust")
source(here::here("code", "make_functions.R"))


###################################################
#                  Parameters                     #
###################################################
# To select which tables are included in the report.
# Also to modify the headers and footers for each table.
tlf <-
  list(
    tab_dm_neg = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the baseline SARS-CoV-2 negative per-protocol cohort",
      table_footer = "This table summarises the per-protocol individuals 
      without a COVID failure event $<$ 7 days post Day 57.",
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_dm_pos = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the baseline SARS-CoV-2 positive per-protocol cohort",
      table_footer = "This table summarises the per-protocol individuals 
      without a COVID failure event $<$ 7 days post Day 57.",
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    case_vacc_neg = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2 negative
      per-protocol cohort (vaccine recipients)",
      table_footer =
        "*Cases are baseline negative per-protocol vaccine recipients with the
        symptomatic infection COVID-19 primary endpoint diagnosed starting 7 days
        after the Day 57 study visit.  Non-cases/Controls are baseline negative
        per-protocol vaccine recipients sampled into the random subcohort with
        no evidence of SARS-CoV-2 infection up to the time of data cut.",
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases*" = 3, "Non-Cases/Control" = 3,
                        "Comparison" = 2),
      header_above2 = c(" "=2,
                        "Baseline SARS-CoV-2 Negative Vaccine Recipients" = 8),
      col1="1cm"
    ),
    
    case_vacc_pos = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2 positive
      per-protocol cohort (vaccine recipients)",
      table_footer = c(
        "*Cases are baseline positive per-protocol vaccine recipients with the
        symptomatic infection COVID-19 primary endpoint diagnosed starting 7
        days after the Day 57 study visit.  Non-cases/Controls are baseline
        negative per-protocol vaccine recipients sampled into the random
        subcohort with no evidence of SARS-CoV-2 infection up to the time
        of data cut."),
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases*" = 3, "Non-Cases/Control" = 3, 
                        "Comparison" = 2),
      header_above2 = c(" "=2,
                        "Baseline SARS-CoV-2 Positive Vaccine Recipients" = 8),
      col1="1cm"
    ),
    
    case_plcb_pos = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2 positive
      per-protocol cohort (placebo recipients)",
      table_footer = c(
        "*Cases are baseline negative per-protocol vaccine recipients with the
        symptomatic infection COVID-19 primary endpoint diagnosed starting 7
        days after the Day 57 study visit.  Non-cases/Controls are baseline
        negative per-protocol vaccine recipients sampled into the random
        subcohort with no evidence of SARS-CoV-2 infection up to the time of
        data cut."),
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2,  "Cases*" = 3, "Non-Cases/Control" = 3,
                        "Comparison" = 2),
      header_above2 = c(" "=2,
                        "Baseline SARS-CoV-2 Positive Placebo Recipients" = 8),
      col1="1cm"
    )
  )


# Depends on the Incoming data

llods <-c(bindN = 20, bindSpike = 20, bindRBD = 20, pseudoneutid50 = 10, 
          pseudoneutid80 = 10, liveneutmn50 = 62.16) 
lloqs <-c(bindN = 34, bindSpike = 34, bindRBD = 34, pseudoneutid50 = 49, 
          pseudoneutid80 = 43, liveneutmn50 = 117.35) 
uloqs <-c(bindN = 19136250, bindSpike = 19136250, bindRBD = 19136250, 
          pseudoneutid50 = Inf, pseudoneutid80 = Inf, liveneutmn50 = 18976.19) 

labels.assays.short <- c(bindN = "Anti N IgG (IU/ml)", 
                         bindSpike = "Anti Spike IgG (IU/ml)", 
                         bindRBD = "Anti RBD IgG (IU/ml)", 
                         pseudoneutid50 = "Pseudovirus-nAb ID50", 
                         pseudoneutid80 = "Pseudovirus-nAb ID80", 
                         liveneutmn50 = "Live virus-nAb MN50")

labels.time <- c(B = "Day 1", Day29 = "Day 29", Day57 = "Day 57", 
                 Delta29overB = "D29 fold-rise over D1", 
                 Delta57overB = "D57 fold-rise over D1", 
                 Delta57over29 = "D57 fold-rise over D29")

assays <- unique(c("bindN"[include_bindN], assays))
labels.assays.short <- labels.assays.short[assays]
labels.time <- labels.time[times]

# 

labels.assays.long <- data.frame (purrr::imap_dfc(labels.assays.short, ~ paste0(labels.assays.short[.y], ": ", labels.time)))
rownames(labels.assays.long) <- names(labels.time)


visits <- names(labels.time)[!grepl("Delta", names(labels.time))]
assays_col <- levels(interaction(visits, assays, sep=""))

labels.assays <- expand.grid(
  time = rownames(labels.assays.long),
  marker = colnames(labels.assays.long),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    label.long = labels.assays.long[time, marker],
    label.short = sapply(labels.assays.short, as.character)[marker],
    Marker = strsplit(as.character(label.long), ": ", fixed = T)[[1]][1],
    Visit = strsplit(as.character(label.long), ": ", fixed = T)[[1]][2],
    colname = paste0(time, marker)
  )

resp.lb <- expand.grid(
  time = visits, marker = assays,
  ind = c("Resp", "FR2", "FR4", "2llod", "4llod"), stringsAsFactors = F
) %>%
  mutate(Ind = case_when(
    ind == "FR2" ~ "% 2-Fold Rise",
    ind == "FR4" ~ "% 4-Fold Rise",
    ind == "Resp" ~ "Responder",
    ind == "2llod" ~ "% Greater than 2xLLOD",
    ind == "4llod" ~ "% Greater than 4xLLOD"
  )) 

labels_all <- full_join(labels.assays, resp.lb, by = c("time", "marker")) %>% 
  mutate(mag_cat = colname, resp_cat = paste0(colname, ind))


###################################################
#                Clean the Data                   #
###################################################

### Table 1. Demographics 
# Output: tab_dm
# Select the covariates to be summarised.
# num_v are columns from ds_long;
# cat_v are rows of `subgroup`


# Read in original data
dat.mock <- read.csv(here::here("../data_clean", data_name))

# The stratified random cohort for immunogenicity
ds_s <- dat.mock %>%
  # dplyr::filter(Perprotocol == 1 & EventTimePrimaryD57 >= 7) %>%
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
    cohort1 = (Perprotocol == 1 & EventTimePrimaryD57 >= 7),
    corrset1 = !is.na(wt))

# SAP Section 9.1.1 Antibody Marker Data
# Step1: Truncation - LLOD, ULOQ
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

ds <- bind_cols(
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

subgrp <- c(
  All = "All participants", 
  Age65C = "Age",
  BMI="BMI",
  HighRiskC = "Risk for Severe Covid-19",
  AgeRiskC = "Age, Risk for Severe Covid-19",
  AgeRisk1 = "Age $<$ 65, Risk for Severe Covid-19",
  AgeRisk2 = "Age $\\geq$ 65, Risk for Severe Covid-19",
  SexC = "Sex", 
  AgeSexC = "Age, sex",
  ethnicityC = "Hispanic or Latino ethnicity", 
  RaceEthC = "Race",
  MinorityC = "Underrepresented minority status",
  AgeMinorC = "Age, Underrepresented minority status"

)

###################################################
#             Generating the Tables               #
###################################################

num_v1 <- c("Age") # Summaries - Mean & Range
num_v2 <- c("BMI") # Summaries - Mean & St.d
cat_v <- c("Age65C", "SexC", "RaceEthC", "ethnicityC", "HighRiskC", "AgeRiskC")

ds_long_ttl <- ds %>%
  dplyr::filter(corrset1) %>% 
  bind_rows(mutate(., Arm="Total")) %>% 
  mutate(AgeRiskC = ifelse(grepl("$\\geq$", AgeRiskC, fixed=T), "Age $\\geq$ 65 ", AgeRiskC),
         Day=57) %>% 
  mutate_all(as.character) %>% 
  pivot_longer(all_of(c(num_v1, num_v2, cat_v)), names_to="subgroup", values_to="subgroup_cat")

if(has29){
  ds <- mutate(ds, 
               cohort2 = (EventTimePrimaryD29 >= 14 & Perprotocol == 1 |
                          EventTimePrimaryD29 >= 7 & EventTimePrimaryD29 <= 13 & Fullvaccine == 1),
               corrset2 = !is.na(wt.2))

  ds_long_ttl <- bind_rows(
    ds_long_ttl, 
    ds %>%
     dplyr::filter(corrset2) %>% 
     bind_rows(mutate(., Arm="Total")) %>% 
     mutate(AgeRiskC = ifelse(grepl("$\\geq$", AgeRiskC, fixed=T), "Age $\\geq$ 65 ", AgeRiskC),
            Day=29) %>% 
     mutate_all(as.character) %>% 
     pivot_longer(all_of(c(num_v1, num_v2, cat_v)), names_to="subgroup", values_to="subgroup_cat"))
}

# Calculate % for categorical covariates
dm_cat <- inner_join(
  ds_long_ttl %>%
    group_by(Day, `Baseline SARS-CoV-2`, Arm, subgroup, subgroup_cat) %>%
    summarise(n = n(), .groups = 'drop'),
  ds_long_ttl %>%
    group_by(Day, `Baseline SARS-CoV-2`, Arm, subgroup) %>%
    summarise(N = n(), .groups = 'drop'),
  by = c("Day", "Baseline SARS-CoV-2", "Arm", "subgroup")
) %>%
  mutate(pct = n / N,
         rslt1 = sprintf("%s (%.1f%%)", n, n / N * 100), 
         rslt2 = sprintf("%s/%s = %.1f%%", n, N, n / N * 100),
         subgroup = ifelse(subgroup_cat == "Communities of Color", 
                           "Race", as.character(subgroup))) %>% 
  dplyr::filter(subgroup %in% cat_v) 

# Calculate mean and range for numeric covariates
dm_num <- ds_long_ttl %>%
  dplyr::filter(subgroup %in% c(num_v1, num_v2)) %>% 
  mutate(subgroup_cat=as.numeric(subgroup_cat)) %>%
  group_by(Day, `Baseline SARS-CoV-2`, Arm, subgroup) %>%
  summarise(
    min = min(subgroup_cat, na.rm = T), 
    max = max(subgroup_cat, na.rm = T),
    mean = mean(subgroup_cat, na.rm = T),
    sd = sd(subgroup_cat, na.rm = T), 
    rslt1 = sprintf("%.1f (%.1f, %.1f)", mean, min, max),
    rslt2 = sprintf("%.1f $\\pm$ %.1f", mean, sd),
    N = n(),
    .groups = 'drop') %>% 
  mutate(subgroup_cat = case_when(subgroup %in% num_v1 ~ "Mean (Range)",
                                  subgroup %in% num_v2 ~ "Mean $\\pm$ SD"),
         subgroup=ifelse(subgroup=="Age", "Age65C", subgroup))

char_lev <- c("Age $<$ 65", "Age $\\geq$ 65", "Mean (Range)","Mean $\\pm$ SD",
              "Female","Male","White Non-Hispanic ","Black or African American",
              "Asian", "American Indian or Alaska Native",
              "Native Hawaiian or Other Pacific Islander", "Multiracial",
              "Other", "Not reported and unknown", "Communities of Color",
              "Hispanic or Latino","Not Hispanic or Latino",
              "Not reported and unknown ","At-risk","Not at-risk",
              "Age $<$ 65 At-risk","Age $<$ 65 Not at-risk", "Age $\\geq$ 65 ")


tab_dm <- bind_rows(dm_cat, dm_num) %>%
  mutate(rslt = case_when(subgroup %in% cat_v ~ rslt1,
                          subgroup %in% num_v1 ~ rslt1,
                          subgroup %in% num_v2 ~ rslt2)) %>%
  dplyr::filter(subgroup_cat %in% char_lev) %>% 
  inner_join(ds_long_ttl %>% 
               distinct(Day, `Baseline SARS-CoV-2`, Arm, Ptid) %>% 
               group_by(Day, `Baseline SARS-CoV-2`, Arm) %>%
               summarise(tot = n()),
             by = c("Day", "Baseline SARS-CoV-2", "Arm")) %>% 
  mutate(Arm = paste0(Arm, "\n(N = ", tot, ")"), subgroup=subgrp[subgroup]) %>%
  pivot_wider(c(Day, `Baseline SARS-CoV-2`, Arm, subgroup, subgroup_cat, rslt),
              names_from = Arm,
              names_sort = T,
              values_from = c(rslt)) %>%
  mutate(Characteristics = factor(subgroup_cat, levels=char_lev), 
         subgroup=factor(subgroup, subgrp)) %>%
  arrange(Day, `Baseline SARS-CoV-2`, subgroup, Characteristics)

tab_dm_pos <- tab_dm %>% 
  dplyr::filter(`Baseline SARS-CoV-2` == "Positive" & Day == 57) %>% 
  select_if(~ !all(is.na(.))) %>% 
  select_at(c("subgroup", "Characteristics", 
              grep("Vaccine" ,names(.), value = T),
              grep("Placebo" ,names(.), value = T),
              grep("Total" ,names(.), value = T)))

tab_dm_neg <- tab_dm %>% 
  dplyr::filter(`Baseline SARS-CoV-2` == "Negative" & Day == 57) %>% 
  select_if(~ !all(is.na(.))) %>% 
  select_at(c("subgroup", "Characteristics", 
              grep("Vaccine" ,names(.), value = T),
              grep("Placebo" ,names(.), value = T),
              grep("Total" ,names(.), value = T)))

if(has29){
  tab_dm_pos.2 <- tab_dm %>% 
    dplyr::filter(`Baseline SARS-CoV-2` == "Positive" & Day == 29) %>% 
    select_if(~ !all(is.na(.))) %>% 
    select_at(c("subgroup", "Characteristics", 
                grep("Vaccine" ,names(.), value = T),
                grep("Placebo" ,names(.), value = T),
                grep("Total" ,names(.), value = T)))
  
  tab_dm_neg.2 <- tab_dm %>% 
    dplyr::filter(`Baseline SARS-CoV-2` == "Negative" & Day == 29) %>% 
    select_if(~ !all(is.na(.))) %>% 
    select_at(c("subgroup", "Characteristics", 
                grep("Vaccine" ,names(.), value = T),
                grep("Placebo" ,names(.), value = T),
                grep("Total" ,names(.), value = T)))
} else {
  tab_dm_pos.2 <- tab_dm_neg.2 <- NULL
}


print("Done with table 1") 

# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021)
# Cases vs Non-cases

corr.design1 <- twophase(list(~Ptid, ~Ptid), 
                        strata=list(NULL, ~ Wstratum),
                        weights=list(NULL, ~ wt),
                        subset= ~ corrset1,
                        method="simple",
                        data=subset(ds, cohort1))

sub.by <- c("Arm", "`Baseline SARS-CoV-2`")
resp.v <- grep("Resp", names(ds), value = T) 
gm.v <- assays_col
subs=c("Case")
comp_i <- c("Cases", "Non-Cases")

rpcnt_case <- get_rr(ds, resp.v, subs, sub.by, design=corr.design1, "wt", "corrset1") %>% 
  mutate(Day=57)
rgm_case <- get_gm(ds, gm.v, subs, sub.by, design=corr.design1, "corrset1") %>% 
  mutate(Day=57)
rgmt_case <- get_rgmt(ds, gm.v, subs, comp_lev=comp_i, sub.by, "Wstratum", "wt", "corrset1") %>% 
  mutate(Day=57)

print("Done with table 2 & 3") 

if(has29){
  corr.design2 <- twophase(list(~Ptid, ~Ptid), 
                           strata=list(NULL, ~ Wstratum),
                           weights=list(NULL, ~ wt.2),
                           subset= ~ corrset2,
                           method="simple",
                           data=subset(ds, cohort2))
  
  rpcnt_case2 <- get_rr(ds, resp.v, subs, sub.by, design=corr.design2, "wt.2", "corrset2")
  rgm_case2 <- get_gm(ds, gm.v, subs, sub.by, design=corr.design2, "corrset2")
  rgmt_case2 <- get_rgmt(ds, gm.v, subs, comp_lev=comp_i, sub.by, "Wstratum", "wt.2", "corrset2")

  rpcnt_case <- bind_rows(rpcnt_case, mutate(rpcnt_case2, Day=29))
  rgm_case <- bind_rows(rgm_case, mutate(rgm_case2, Day=29))
  rgmt_case <- bind_rows(rgmt_case, mutate(rgmt_case2, Day=29))
  
  print("Done with table 2b & 3b") 
}

rrdiff_case <- rpcnt_case %>% 
  dplyr::filter(subgroup %in% subs & grepl("Resp",resp_cat)) %>% 
  mutate(groupn = 2-match(Group, comp_i)%%2) %>%
  pivot_wider(id_cols = c(Day, subgroup, `Baseline SARS-CoV-2`, Arm, Visit, Marker, Ind),
              names_from = groupn, values_from = c(response, ci_l, ci_u), names_sep = "") %>% 
  mutate(Estimate = response1-response2,
         ci_l = Estimate-sqrt((response1-ci_l1)^2+(response2-ci_u2)^2),
         ci_u = Estimate+sqrt((response1-ci_u1)^2+(response2-ci_l2)^2),
         rrdiff = sprintf("%s%%\n(%s%%, %s%%)", round(Estimate * 100, 1), 
                          round(ci_l*100, 1), round(ci_u*100, 1))) 

print("Done with table6") 

tab_case <- inner_join(rpcnt_case, rgm_case,
                       by = c("Day", "Group", "Arm", "Baseline SARS-CoV-2", 
                              "N", "Marker", "Visit")) %>%
  pivot_wider(id_cols = c(Day, Arm, `Baseline SARS-CoV-2`, Marker, Visit),
              names_from = Group, 
              values_from = c(N, rslt, `GMT/GMC`)) %>% 
  inner_join(rrdiff_case, by = c("Day", "Arm", "Baseline SARS-CoV-2", "Marker", "Visit")) %>% 
  inner_join(rgmt_case, by = c("Day", "Arm", "Baseline SARS-CoV-2", "Marker", "Visit")) %>% 
  filter((Day==57 & Visit=="Day 57")|(Day==29 & Visit=="Day 29")) %>% 
  select(Arm, `Baseline SARS-CoV-2`, Visit, Marker, `N_Cases`, `rslt_Cases`, 
         `GMT/GMC_Cases`, `N_Non-Cases`, `rslt_Non-Cases`, `GMT/GMC_Non-Cases`,  
         rrdiff, `Ratios of GMT/GMC`) %>% 
  arrange(Arm, `Baseline SARS-CoV-2`, Visit) 
  
case_vacc_neg <- tab_case %>% 
  dplyr::filter(Arm == "Vaccine" & `Baseline SARS-CoV-2` == "Negative") %>% 
  select(-c(Arm, `Baseline SARS-CoV-2`))

case_vacc_pos <- tab_case %>% 
  dplyr::filter(Arm == "Vaccine" & `Baseline SARS-CoV-2` == "Positive") %>% 
  select(-c(Arm, `Baseline SARS-CoV-2`))

case_plcb_pos <- tab_case %>% 
  dplyr::filter(Arm == "Placebo" & `Baseline SARS-CoV-2` == "Positive") %>% 
  select(-c(Arm, `Baseline SARS-CoV-2`))

print("Done with all tables") 

save(tlf, tab_dm_neg, tab_dm_pos, tab_dm_neg.2, tab_dm_pos.2, 
     case_vacc_neg, case_vacc_pos, case_plcb_pos,
     file = here::here("output", "Tables.Rdata"))


