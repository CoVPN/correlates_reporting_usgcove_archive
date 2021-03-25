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
      table_footer = "This table summarises the random subcohort, which
      was randomly sampled from the per-protocol individuals without a COVID failure
      event < 7 days post Day 57. The sampling was stratified by the key baseline 
      covariates: assigned treatment arm, baseline SARS-CoV-2 status 
      (defined by serostatus and possibly also NAAT and/or RNA PCR testing), 
      any additional important demographic factors such as the randomization strata 
      (e.g., defined by age and/or co-morbidities).",
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_dm_pos = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the baseline SARS-CoV-2 positive per-protocol cohort",
      table_footer = "This table summarises the random subcohort, which
      was randomly sampled from the per-protocol individuals without a COVID failure
      event < 7 days post Day 57. The sampling was stratified by the key baseline 
      covariates: assigned treatment arm, baseline SARS-CoV-2 status 
      (defined by serostatus and possibly also NAAT and/or RNA PCR testing), 
      any additional important demographic factors such as the randomization strata 
      (e.g., defined by age and/or co-morbidities).",
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

ds <- ds2

resp_v <- setdiff(grep("Resp", names(ds), value = T),
                  grep("liveneut", names(ds), value = T)
)

ds_resp_l <- pivot_longer(ds,
                          cols = all_of(resp_v), 
                          names_to = "resp_cat",
                          values_to = "response") %>%
  inner_join(labels_all, by = "resp_cat")

mag_v <- c(assays_col, grep("Delta", names(ds), value = T))

ds_mag_l <- pivot_longer(ds,
                         cols = all_of(mag_v),
                         names_to = "mag_cat", values_to = "mag") %>%
  inner_join(labels_all %>% 
               distinct(mag_cat, time, marker, Visit, Marker, label.short), 
             by = "mag_cat")


###################################################
#             Generating the Tables               #
###################################################

num_v1 <- c("Age") # Summaries - Mean & Range
num_v2 <- c("BMI") # Summaries - Mean & St.d
cat_v <- c("Age", "Risk for Severe Covid-19", "Sex", "Race", 
           "Hispanic or Latino ethnicity", "Risk for Severe Covid-19")

# Stack a Arm = "Total" to the original data
ds_long_ttl <- bind_rows(
  ds_long %>% mutate(Arm = "Total"),
  ds_long
)

# Calculate % for categorical covariates
dm_cat <- inner_join(
  ds_long_ttl %>%
    group_by(`Baseline SARS-CoV-2`, Arm, subgroup, subgroup_cat) %>%
    summarise(n = n(), .groups = 'drop'),
  ds_long_ttl %>%
    group_by(`Baseline SARS-CoV-2`, Arm, subgroup) %>%
    summarise(N = n(), .groups = 'drop'),
  by = c("Baseline SARS-CoV-2", "Arm", "subgroup")
) %>%
  mutate(pct = n / N,
         rslt1 = sprintf("%s (%.1f%%)", n, n / N * 100), 
         rslt2 = sprintf("%s/%s = %.1f%%", n, N, n / N * 100),
         subgroup = ifelse(subgroup_cat == "Communities of Color", 
                           "Race", as.character(subgroup))) %>% 
  dplyr::filter(subgroup %in% cat_v) %>% 
  arrange(`Baseline SARS-CoV-2`, Arm, subgroup)

# Calculate mean and range for numeric covariates
dm_num <- ds_long_ttl %>%
  distinct_at(all_of(c("Baseline SARS-CoV-2","Arm","Ptid",num_v1,num_v2))) %>%
  pivot_longer(cols = all_of(c(num_v1, num_v2)), 
               names_to = "subgroup", 
               values_to = "subgroup_cat") %>%
  group_by(`Baseline SARS-CoV-2`, Arm, subgroup) %>%
  summarise(
    min = min(subgroup_cat, na.rm = T), 
    max = max(subgroup_cat, na.rm = T),
    mean = mean(subgroup_cat, na.rm = T),
    sd = sd(subgroup_cat, na.rm = T), 
    rslt1 = sprintf("%.1f (%.1f, %.1f)", mean, min, max),
    rslt2 = sprintf("%.1f $\\pm$ %.1f", mean, sd),
    N = n(),
    .groups = 'drop'
  ) %>% 
  mutate(subgroup_cat = case_when(subgroup %in% num_v1 ~ "Mean (Range)",
                                  subgroup %in% num_v2 ~ "Mean $\\pm$ SD"))

tab_dm <- full_join(
  dm_cat, 
  dm_num,
  by = c("Baseline SARS-CoV-2", "Arm", "subgroup", 
         "subgroup_cat", "N", "rslt1", "rslt2")
) %>%
  mutate(rslt = case_when(subgroup %in% cat_v ~ rslt1,
                          subgroup %in% num_v1 ~ rslt1,
                          subgroup %in% num_v2 ~ rslt2),
         subgroup = factor(subgroup, 
                           levels = c("Age", "BMI", "Sex", "Race", 
                                      "Hispanic or Latino ethnicity", 
                                      "Risk for Severe Covid-19"))) %>%
  inner_join(ds_long_ttl %>% 
               distinct(`Baseline SARS-CoV-2`, Arm, Ptid) %>% 
               group_by(`Baseline SARS-CoV-2`, Arm) %>%
               summarise(tot = n()),
             by = c("Baseline SARS-CoV-2", "Arm")) %>% 
  mutate(Arm = paste0(Arm, "\n(N = ", tot, ")")) %>%
  pivot_wider(c(`Baseline SARS-CoV-2`, Arm, subgroup, subgroup_cat, rslt),
              names_from = Arm,
              names_sort = T,
              values_from = c(rslt)) %>%
  mutate(Characteristics = factor(subgroup_cat,
                                  levels=c("Age $\\geq$ 65","Age $<$ 65",
                                           "Mean (Range)","Mean $\\pm$ SD",
                                           "Female","Male",
                                           "White Non-Hispanic ",
                                           "Black or African American",
                                           "Asian",
                                           "American Indian or Alaska Native",
                                           "Native Hawaiian or Other Pacific Islander",
                                           "Multiracial",
                                           "Other",
                                           "Not reported and unknown",
                                           "Communities of Color",
                                           "Hispanic or Latino","Not Hispanic or Latino",
                                           "Not reported and unknown ",
                                           "At-risk","Not at-risk"))) %>% 
  arrange(`Baseline SARS-CoV-2`, subgroup, Characteristics)

tab_dm_pos <- tab_dm %>% 
  dplyr::filter(`Baseline SARS-CoV-2` == "Positive") %>% 
  select_if(~ !all(is.na(.))) %>% 
  select_at(c("subgroup", "Characteristics", 
                         grep("Vaccine" ,names(.), value = T),
                         grep("Placebo" ,names(.), value = T),
                         grep("Total" ,names(.), value = T)))

tab_dm_neg <- tab_dm %>% 
  dplyr::filter(`Baseline SARS-CoV-2` == "Negative") %>% 
  select_if(~ !all(is.na(.))) %>% 
  select_at(c("subgroup", "Characteristics", 
              grep("Vaccine" ,names(.), value = T),
              grep("Placebo" ,names(.), value = T),
              grep("Total" ,names(.), value = T)))


print("Done with table 1") 

# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021)
# Cases vs Non-cases

ds_resp_case <- ds_resp_l %>% 
  dplyr::filter(subgroup == "All participants" & grepl("Resp", resp_cat))

resp_v <- unique(ds_resp_case$resp_cat)

dssvy_case <- svydesign(ids = ~ Ptid, 
                        strata = ~ Wstratum,
                        weights = ~ wt,
                        data = ds_resp_case,
                        nest = F)


sub_grp_col <- c("subgroup", "Arm", "`Baseline SARS-CoV-2`", 
                 "Case", "resp_cat", "Visit", "Marker")
sub_grp <- as.formula(paste0("~", paste(sub_grp_col, collapse = "+")))

rpcnt_case <- svyby(~response, by=sub_grp, dssvy_case, svyciprop, vartype = "ci")

tab_rr_case <- rpcnt_case %>%
  inner_join(
    ds_resp_case %>%
      mutate(rspndr = response*wt.subcohort) %>% 
      group_by(across(all_of(gsub("`","",sub_grp_col)))) %>%
      summarise(N = n(), Nw = sum(wt.subcohort), rspndr = sum(rspndr),
                .groups = 'drop'),
    by = gsub("`","",sub_grp_col)
  ) %>%
  mutate(Responder = case_when(
    is.na(ci_l)|is.na(ci_u) ~ 
      sprintf("%s/%s = %.1f%%", round(rspndr,1), round(Nw,1), response*100),
    TRUE ~ 
      sprintf("%s/%s = %.1f%%\n(%.1f%%, %.1f%%)", 
              round(rspndr, 1), round(Nw, 1), response*100, ci_l*100, ci_u*100))
  ) %>% 
  select(subgroup, Arm, `Baseline SARS-CoV-2`, Case, Visit, N, Marker, Responder)

#########

# 8b Responder rate differences between cases vs non-cases & 95% CI of 
# Titers or Concentrations
comp_v <- c("Cases", "Non-Cases")

rrdiff_case <- rpcnt_case %>% 
  mutate(Case = match(as.character(Case), comp_v), comp = paste(comp_v, collapse = " vs ")) %>% 
  pivot_wider(names_from = Case, values_from = c(response, ci_l, ci_u), names_sep = "") %>% 
  mutate(Estimate = response1-response2,
         ci_l = Estimate-sqrt((response1-ci_l1)^2+(response2-ci_u2)^2),
         ci_u = Estimate+sqrt((response1-ci_u1)^2+(response2-ci_l2)^2)) %>% 
  select(-c(response1, response2, ci_l1, ci_l2, ci_u1, ci_u2))
  

tab_rrdiff_case <- rrdiff_case %>% 
  mutate(rslt = case_when(
    !complete.cases(ci_l, ci_u) ~ sprintf("%s%%", round(Estimate * 100, 1)),
    complete.cases(ci_l, ci_u) ~ sprintf("%s%%\n(%s%%, %s%%)", round(Estimate * 100, 1), 
                                         round(ci_l*100, 1), round(ci_u*100, 1)))
  )

print("Done with table9") 

#########
gm_v <- assays_col
ds_mag_case <- ds_mag_l %>% 
  dplyr::filter(subgroup == "All participants" & mag_cat %in% gm_v & !is.na(wt))

sub_grp_col <- c("subgroup", "Arm", "`Baseline SARS-CoV-2`", "Case", "mag_cat")
sub_grp <- as.formula(paste0("~", paste(sub_grp_col, collapse = "+")))

dssvy_case <- svydesign(ids = ~ Ptid, 
                        strata = ~ Wstratum,
                        weights = ~ wt,
                        data = ds_mag_case,
                        nest = F)

gm_f <- paste0("~", paste(gm_v, collapse = " + "))
rgm_case <- svyby(~ mag, sub_grp, dssvy_case, svymean, vartype = "ci")

tab_gm_case <- rgm_case %>%
  inner_join(ds_mag_case %>%
               group_by(subgroup, mag_cat, Case, Visit, Marker, Arm, `Baseline SARS-CoV-2`) %>%
               summarise(N = n(), .groups = 'drop'),
             by = c("subgroup", "Case", "Arm", "Baseline SARS-CoV-2", "mag_cat")
  ) %>%
  mutate(`GMT/GMC` = sprintf("%.0f\n(%.0f, %.0f)", 10^mag, 10^ci_l, 10^ci_u)
  ) %>% 
  select(subgroup, Arm, `Baseline SARS-CoV-2`, Case, Visit, N, Marker, `GMT/GMC`)

###
ds_mag_l_rgmt_case <- ds_mag_l %>%
  dplyr::filter(subgroup == "All participants" & !is.na(wt) & 
                  !grepl("Delta", mag_cat))

comp_v <- "Case"
comp_lev <- c("Cases", "Non-Cases")
f_v <- as.formula(sprintf("mag ~ %s", comp_v))
sub_grp_col <- c("subgroup", "Arm", "Baseline SARS-CoV-2", "Visit", "Marker")

rgmt_case <- ds_mag_l_rgmt_case %>%
  group_split(across(all_of(sub_grp_col))) %>%
  map_dfr(get_rgmt, comp_v=comp_v, comp_lev=comp_lev, f_v = f_v, 
          sub_grp_col = sub_grp_col, stratum = "Wstratum", weights = "wt") 

tab_rgmt_case <- rgmt_case %>%
  mutate(`Ratios of GMT/GMC` = sprintf("%.2f\n(%.2f, %.2f)", 
                                       10^Estimate, 10^ci_l, 10^ci_u)) %>% 
  select(Arm, `Baseline SARS-CoV-2`, Marker, Visit, subgroup, comp, `Ratios of GMT/GMC`) %>% 
  arrange(subgroup, Arm, `Baseline SARS-CoV-2`, Visit, Marker)

tab_key_case <- full_join(tab_rr_case, 
                          tab_gm_case,
                          by = c("subgroup", "Arm", "Baseline SARS-CoV-2", 
                                 "Case", "N", "Marker", "Visit")) %>%
  pivot_wider(names_from = Case, 
              values_from = c(N, Responder, `GMT/GMC`)) %>% 
  # Join with RR difference (tab_rr)
  full_join(tab_rrdiff_case, 
            by = c("Arm", "subgroup", "Baseline SARS-CoV-2", "Visit", "Marker")
  ) %>% 
  # Join with GMT/GMC Ratios (tab_rgmt) 
  full_join(tab_rgmt_case, 
            by = c("Arm", "Baseline SARS-CoV-2", "Visit", "Marker", "subgroup")
  ) %>% 
  select(Arm, `Baseline SARS-CoV-2`, Visit, Marker, `N_Cases`, `Responder_Cases`, 
         `GMT/GMC_Cases`, `N_Non-Cases`, `Responder_Non-Cases`, `GMT/GMC_Non-Cases`,  
         rslt, `Ratios of GMT/GMC`) %>% 
  dplyr::filter(Visit != "Day 1") %>% 
  arrange(Arm, `Baseline SARS-CoV-2`) 

case_vacc_neg <- filter(tab_key_case, Arm == "Vaccine", `Baseline SARS-CoV-2` == "Negative") %>% 
  select(-c(Arm, `Baseline SARS-CoV-2`))
case_vacc_pos <- filter(tab_key_case, Arm == "Vaccine", `Baseline SARS-CoV-2` == "Positive") %>% 
  select(-c(Arm, `Baseline SARS-CoV-2`))
case_plcb_pos <- filter(tab_key_case, Arm == "Placebo", `Baseline SARS-CoV-2` == "Positive") %>% 
  select(-c(Arm, `Baseline SARS-CoV-2`))

print("Done with table10-12") 


save(tlf, tab_dm_neg, tab_dm_pos, case_vacc_neg, case_vacc_pos, case_plcb_pos, 
     file = here::here("output", "Tables.Rdata"))


