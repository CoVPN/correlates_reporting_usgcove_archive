##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here("..")) #
source(here::here("..", "_common.R")) #
##################################################

# Immunogenicity Tables

# Reload clean_data
base::load(here::here("data_clean", "ds_all.Rdata"))
source(here::here("code", "make_functions.R"))
library(survey)
library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
# For stratum with 1 ppt
options(survey.lonely.psu="adjust")
###################################################
#             Generating the Tables               #
###################################################

### Table 1. Demographics 
# Output: tab_dm
# Select the covariates to be summarised.
# num_v are columns from ds_long;
# cat_v are rows of `subgroup`


num_v1 <- c("Age") # Summaries - Mean & Range
num_v2 <- c("BMI") # Summaries - Mean & St.d
cat_v <- c("Age65C", "SexC", "RaceEthC", "ethnicityC", "HighRiskC", "AgeRiskC")

ds_long_ttl <- ds %>%
  dplyr::filter(randomset) %>% 
  bind_rows(mutate(., Arm="Total")) %>% 
  mutate(AgeRiskC = ifelse(grepl("$\\geq$", AgeRiskC, fixed=T), "Age $\\geq$ 65 ", AgeRiskC)) %>% 
  mutate_all(as.character) %>% 
  pivot_longer(all_of(c(num_v1, num_v2, cat_v)), names_to="subgroup", values_to="subgroup_cat")

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
  dplyr::filter(subgroup %in% c(num_v1, num_v2)) %>% 
  mutate(subgroup_cat=as.numeric(subgroup_cat)) %>%
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

char_lev <- c("Age $<$ 65", "Age $\\geq$ 65", "Mean (Range)","Mean $\\pm$ SD",
              "Female","Male","White Non-Hispanic ","Black or African American",
              "Asian", "American Indian or Alaska Native",
              "Native Hawaiian or Other Pacific Islander", "Multiracial",
              "Other", "Not reported and unknown", "Communities of Color",
              "Hispanic or Latino","Not Hispanic or Latino",
              "Not reported and unknown ","At-risk","Not at-risk",
              "Age $<$ 65 At-risk","Age $<$ 65 Not at-risk", "Age $\\geq$ 65 ")

tab_dm <- full_join(dm_cat, dm_num,
  by = c("Baseline SARS-CoV-2", "Arm", "subgroup", 
         "subgroup_cat", "N", "rslt1", "rslt2")
) %>%
  mutate(rslt = case_when(subgroup %in% cat_v ~ rslt1,
                          subgroup %in% num_v1 ~ rslt1,
                          subgroup %in% num_v2 ~ rslt2)) %>%
  dplyr::filter(subgroup_cat %in% char_lev) %>% 
  inner_join(ds_long_ttl %>% 
               distinct(`Baseline SARS-CoV-2`, Arm, Ptid) %>% 
               group_by(`Baseline SARS-CoV-2`, Arm) %>%
               summarise(tot = n()),
             by = c("Baseline SARS-CoV-2", "Arm")) %>% 
  mutate(Arm = paste0(Arm, "\n(N = ", tot, ")"), subgroup=subgrp[subgroup]) %>%
  pivot_wider(c(`Baseline SARS-CoV-2`, Arm, subgroup, subgroup_cat, rslt),
              names_from = Arm, 
              names_sort = T,
              values_from = c(rslt)) %>%
  mutate(Characteristics = factor(subgroup_cat, levels=char_lev)) %>%
  arrange(`Baseline SARS-CoV-2`, Characteristics)

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

### Table 2. Responder Rates & Proportions of Magnitudes >= 2FR, 4FR
# For each binding antibody marker, the estimated percentage of participants
# defined as responders, and with concentrations >= 2x LLOD or >=
# 4 x LLOD, will be provided with the corresponding 95% CIs
# 
# Output: tab_bind


# Variables used for stratification in the tables
# subgroup: SAP Table 6: Baseline Subgroups
# Arm and Baseline: Assigned treatment Arms * Baseline SARS-CoV-2-19 Status
# Group: Category in each subgroup

immuno.design <- twophase(list(~Ptid, ~Ptid), 
                          strata=list(NULL, ~tps.stratum),
                          weights=list(NULL, ~wt.subcohort),
                          subset=~randomset,
                          method="simple",
                          data=ds)

sub.by <- c("Arm", "`Baseline SARS-CoV-2`")
resp.v <- grep("Resp|2llod|4llod|FR2|FR4", names(ds), value = T) 
subs=c("All", "Age65C", "HighRiskC", "AgeRiskC", "AgeRisk1", "AgeRisk2", "SexC",
       "AgeSexC", "ethnicityC", "RaceEthC", "MinorityC", "AgeMinorC")
rpcnt <- get_rr(ds, resp.v, subs, sub.by, immuno.design, "wt.subcohort", "randomset")

tab_rr <- rpcnt %>% 
  dplyr::filter(!subgroup %in% c("AgeRisk1", "AgeRisk2") & Visit != "Day 1" & Group %in% names(grplev)) %>% 
  mutate(subgroup=factor(subgrp[subgroup], levels=subgrp), Group=factor(grplev[Group], levels=grplev)) %>% 
  pivot_wider(
    id_cols = c(subgroup, Group, Arm, `Baseline SARS-CoV-2`, Marker, Visit, N),
    names_from = Ind, values_from = rslt) %>% 
  arrange(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker)

if(any(grepl("bind", assays))){
  tab_bind <- tab_rr %>% 
    dplyr::filter(Marker %in% labels_all$Marker[grep("bind", labels_all$marker)]) %>% 
    select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, N, Responder, 
           `% Greater than 2xLLOD`,`% Greater than 4xLLOD`)
}else{
  tab_bind <- NULL
}

print("Done with table2") 

# Table 3 & 4. For the ID50 pseudo-virus & MN50 WT live virus neutralization 
# antibody marker, the estimated percentage of participants defined as 
# responders, participants with % 2-Fold Rise (2FR), and participants with 4-fold 
# rise (4FR) will be provided with the corresponding 95% CIs 
# 
# Output: tab_pseudo & tab_wt

if("pseudoneutid50" %in% assays){
tab_pseudo <- tab_rr %>% 
  dplyr::filter(Marker  %in% labels_all$Marker[grep("pseudoneutid50", labels_all$marker)]) %>% 
  select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, N, Responder,
         `% 2-Fold Rise`, `% 4-Fold Rise`)
}else{
  tab_pseudo <- NULL
}
if("liveneutmn50" %in% assays){
tab_wt <- tab_rr %>% 
  dplyr::filter(Marker  %in% labels_all$Marker[grep("liveneutmn50", labels_all$marker)]) %>% 
  select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, N, Responder,
         `% 2-Fold Rise`, `% 4-Fold Rise`)
}else{
  tab_wt <- NULL
}
print("Done with table3 & 4") 

# Table 5. Geometric mean titers (GMTs) and geometric mean concentrations (GMCs)
# will be summarized along with their 95% CIs using the t-distribution
# approximation of log-transformed concentrations/titers (for each of the 5
# Spike-targeted marker types including pseudovirus-nAb ID50 and ID80
# and WT live virus-nAb MN50, as well as for binding Ab to N).
# 
# Output: tab_gm
gm.v <- c(assays_col, grep("Delta", names(ds), value = T))
rgm <- get_gm(ds, gm.v, subs, sub.by, immuno.design, "randomset")

tab_gm <- rgm %>% 
  dplyr::filter(!subgroup %in% c("AgeRisk1", "AgeRisk2") & !grepl("Delta", mag_cat) & Group %in% names(grplev)) %>% 
  mutate(subgroup=factor(subgrp[subgroup], levels=subgrp), Group=factor(grplev[Group], levels=grplev)) %>% 
  arrange(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker) %>% 
  select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, N, `GMT/GMC`) 

print("Done with table5") 

### Table 6. GMTRs/GMCRs will be summarized with 95% CI (t-distribution 
# approximation) for any post-baseline values compared to baseline, and
# post-Day 57 values compared to Day 57
# 
# Output: tab_gmr

gmr_gm <- inner_join(
  tab_gm %>% 
    dplyr::filter(Visit == "Day 1") %>% 
    select(-Visit) %>% 
    rename(`Baseline\nGMT/GMC` = `GMT/GMC`),
  tab_gm %>% 
    dplyr::filter(Visit != "Day 1") %>% 
    rename(`Post Baseline\nGMT/GMC` = `GMT/GMC`),
  by = c("subgroup", "Arm", "Group", "Baseline SARS-CoV-2", "N", "Marker")) %>% 
  dplyr::filter(!subgroup %in% c("AgeRisk1", "AgeRisk2")) %>% 
  mutate(Visit = paste0(gsub("ay ", "", Visit), " fold-rise over D1"))

tab_gmr <- rgm %>% 
  dplyr::filter(grepl("overB", mag_cat)) %>%
  mutate(`GMTR/GMCR`=sprintf("%.2f\n(%.2f, %.2f)", 10^mag, 10^ci_l, 10^ci_u),
         subgroup=factor(subgrp[subgroup], levels=subgrp), Group=factor(grplev[Group], levels=grplev)) %>% 
  select(subgroup, Group, Arm, `Baseline SARS-CoV-2`, Visit, N, Marker, `GMTR/GMCR`) %>% 
  inner_join(
    gmr_gm,
    c("subgroup", "Group", "Arm", "Baseline SARS-CoV-2", "Visit", "N", "Marker")) %>% 
  arrange(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker) %>% 
  select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, N,  
         `Baseline\nGMT/GMC`, `Post Baseline\nGMT/GMC`, `GMTR/GMCR`) 

print("Done with table6") 

### Table 7. The ratios of GMTs/GMCs will be estimated between groups with the
# two-sided 95% CIs calculated using t-distribution approximation of 
# log-transformed titers/concentrations
# Output: tab_rgmt
# 
# Ratios of GMT/GMC between subgroups among vacinees

comp_lev <- c("Age $\\geq$ 65", "Age $<$ 65",
              "At-risk", "Not at-risk",
              "Age $<$ 65 At-risk", "Age $<$ 65 Not at-risk",
              "Age $\\geq$ 65 At-risk", "Age $\\geq$ 65 Not at-risk",
              "Male", "Female",
              "Hispanic or Latino", "Not Hispanic or Latino",
              "Communities of Color", "White Non-Hispanic")

groups <- c("Age65C", "HighRiskC", "AgeRisk1", "AgeRisk2",
            "SexC", "ethnicityC", "MinorityC")
mag_groups <- assays_col

rgmt <- get_rgmt(ds, mag_groups, groups, comp_lev=comp_lev, 
                 sub.by, "tps.stratum", "wt.subcohort", "randomset")

rgmt_gm <- rgm %>% 
  dplyr::filter(!grepl("Delta", mag_cat) & Group %in% names(grplev)) %>% 
  mutate(subgroup=factor(subgrp[subgroup], levels=subgrp), Group=factor(grplev[Group], levels=grplev)) %>% 
  dplyr::filter(subgroup %in% subgrp[groups]) %>% 
  mutate(groupn = 2-match(Group, comp_lev)%%2) %>% 
  pivot_wider(id_cols = c(subgroup, Arm, `Baseline SARS-CoV-2`, Visit, Marker),
              names_from = groupn, values_from = `GMT/GMC`, 
              names_prefix = "Group")

tab_rgmt <- inner_join(rgmt_gm, 
                       rgmt %>% mutate(subgroup=factor(subgrp[subgroup], levels=subgrp)),  
                       c("Baseline SARS-CoV-2", "Arm", "subgroup", "Marker", "Visit")) %>% 
  rename(`Group 1 vs 2` = comp, 
         `Group 1 GMT/GMC` = `Group1`, 
         `Group 2 GMT/GMC` = `Group2`) %>% 
  select(`Group 1 vs 2`, subgroup, Visit, Arm, `Baseline SARS-CoV-2`, Marker, 
         `Group 1 GMT/GMC`, `Group 2 GMT/GMC`, `Ratios of GMT/GMC`) %>% 
  arrange(subgroup, Visit, Arm, `Baseline SARS-CoV-2`, Marker)

print("Done with table7") 

### Table 8. The differences in the responder rates, 2FRs, 4FRs between groups 
# will be computed along with the two-sided 95% CIs by the Wilson-Score
# method without continuity correction (Newcombe, 1998).
# Output: tab_rrdiff

tab_rrdiff <- bind_rows(rpcnt %>% 
                          dplyr::filter(subgroup=="All") %>% 
                          mutate(Group=Arm, Arm="-"),
                        rpcnt %>% 
                          dplyr::filter(subgroup=="All") %>% 
                          mutate(Group=`Baseline SARS-CoV-2`, `Baseline SARS-CoV-2`="-"),
                        rpcnt)%>% 
  dplyr::filter(subgroup %in% c(groups, "All") & grepl("Resp|FR2|FR4",resp_cat)) %>% 
  mutate(groupn = 2-match(Group, c(comp_lev, "Vaccine", "Placebo", "Positive", "Negative"))%%2) %>% 
  pivot_wider(id_cols = c(subgroup, Arm, `Baseline SARS-CoV-2`, Marker, Visit, Ind),
              names_from = groupn, values_from = c(response, ci_l, ci_u), names_sep = "") %>% 
  full_join(distinct(rgmt, subgroup, comp), by = "subgroup") %>% 
  mutate(Comparison = case_when(Arm=="-"~"Vaccine vs Placebo",
                          `Baseline SARS-CoV-2`=="-"~"Positive vs Negative",
                          TRUE~comp),
         subgroup = factor(case_when(Arm=="-" ~ "Arm",
                           `Baseline SARS-CoV-2`=="-" ~ "Baseline SARS-CoV-2",
                            TRUE~subgrp[subgroup]), levels=c("Arm", "Baseline SARS-CoV-2", subgrp)),
         Estimate = response1-response2,
         ci_l = Estimate-sqrt((response1-ci_l1)^2+(response2-ci_u2)^2),
         ci_u = Estimate+sqrt((response1-ci_u1)^2+(response2-ci_l2)^2),
         rslt = sprintf("%s%%\n(%s%%, %s%%)", 
                        round(Estimate*100,1), round(ci_l*100,1), round(ci_u*100,1))) %>%
  dplyr::filter(!is.na(Comparison)) %>%
  select(Comparison, subgroup, `Baseline SARS-CoV-2`, Arm, Visit, Marker, Ind, rslt) %>%
  pivot_wider(names_from = Ind, values_from = rslt) %>%
  arrange(subgroup, Visit, `Baseline SARS-CoV-2`, Marker, Comparison) 

print("Done with table8") 

# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021) 
# Vaccine vs Placebo within Baseline status
groups <- c("Arm")
comp_i <- c("Vaccine", "Placebo")
mag_groups <- assays_col

rgmt_Rx <- get_rgmt(ds, mag_groups, groups, comp_lev=comp_i, sub.by="`Baseline SARS-CoV-2`", 
                    "tps.stratum", "wt.subcohort", "randomset") %>% 
  mutate(subgroup=factor(subgrp[subgroup], levels=subgrp))

rrdiff_Rx <- rpcnt %>% 
  dplyr::filter(subgroup=="All" & grepl("Resp|FR2|FR4",resp_cat)) %>% 
  mutate(groupn = 2-match(Arm, comp_i)%%2) %>% 
  pivot_wider(id_cols = c(subgroup, `Baseline SARS-CoV-2`, Visit, Marker, Ind),
              names_from = groupn, values_from = c(response, ci_l, ci_u), names_sep = "") %>% 
  mutate(Estimate = response1-response2,
         ci_l = Estimate-sqrt((response1-ci_l1)^2+(response2-ci_u2)^2),
         ci_u = Estimate+sqrt((response1-ci_u1)^2+(response2-ci_l2)^2),
         rrdiff = sprintf("%s%%\n(%s%%, %s%%)", round(Estimate * 100, 1), 
                          round(ci_l*100, 1), round(ci_u*100, 1)),
         subgroup=factor(subgrp[subgroup], levels=subgrp)) 

tab_Rx <- full_join(tab_rr, tab_gm,
                    c("subgroup", "Arm", "Baseline SARS-CoV-2", 
                      "Group", "Visit", "N", "Marker")) %>%
  dplyr::filter(as.character(subgroup) == "All participants") %>%
  pivot_wider(id_cols = c(subgroup, Group, `Baseline SARS-CoV-2`, Marker, Visit),  
              names_from = Arm, 
              values_from = c(N, Responder, `GMT/GMC`)) %>% 
  inner_join(rrdiff_Rx, 
             by = c("subgroup", "Baseline SARS-CoV-2", "Visit", "Marker")) %>% 
  inner_join(rgmt_Rx, 
             by = c("Baseline SARS-CoV-2", "Visit", "Marker")) %>% 
  select(`Baseline SARS-CoV-2`, Visit, Marker, `N_Vaccine`, `Responder_Vaccine`, 
         `GMT/GMC_Vaccine`, `N_Placebo`, `Responder_Placebo`, `GMT/GMC_Placebo`, 
         rrdiff, `Ratios of GMT/GMC`) %>% 
  arrange(`Baseline SARS-CoV-2`, Visit, Marker)

tab_neg <- dplyr::filter(tab_Rx, `Baseline SARS-CoV-2` == "Negative") %>% 
  select(-c(`Baseline SARS-CoV-2`))
tab_pos <- dplyr::filter(tab_Rx, `Baseline SARS-CoV-2` == "Positive") %>% 
  select(-c(`Baseline SARS-CoV-2`))

print("Done with table13-14") 

###################################################

groups <- "`Baseline SARS-CoV-2`"
comp_i <- c("Positive", "Negative")
resp_v <- grep("Resp", names(ds), value=T)

rrdiff_bl <- rpcnt %>% 
  dplyr::filter(subgroup == "All" & grepl("Resp", resp_cat)) %>% 
  mutate(groupn = match(as.character(`Baseline SARS-CoV-2`), comp_i), 
         comp = paste(comp_i, collapse = " vs ")) %>% 
  pivot_wider(id_cols = c(subgroup, Arm, Visit, Marker),
              names_from = groupn, 
              values_from = c(response, ci_l, ci_u), names_sep = "") %>% 
  mutate(Estimate = response1-response2,
         ci_l = Estimate-sqrt((response1-ci_l1)^2+(response2-ci_u2)^2),
         ci_u = Estimate+sqrt((response1-ci_u1)^2+(response2-ci_l2)^2),
         rrdiff = sprintf("%s%%\n(%s%%, %s%%)", round(Estimate * 100, 1), 
                          round(ci_l*100, 1), round(ci_u*100, 1)),
         subgroup=factor(subgrp[subgroup], levels=subgrp))  

rgmt_bl <- get_rgmt(ds, mag_groups, groups, comp_lev = comp_i, sub.by="Arm",
                    "tps.stratum", "wt.subcohort", "randomset") %>% 
  mutate(subgroup=factor(subgrp[subgroup], levels=subgrp))


tab_bl <- full_join(tab_rr, tab_gm, 
                    by = c("Arm", "subgroup", "Group",
                           "Baseline SARS-CoV-2", "Visit", "Marker", "N")) %>%
  dplyr::filter(as.character(subgroup) == "All participants" & Visit!="Day 1") %>%
  pivot_wider(id_cols = c(Arm, subgroup, Group, Marker, Visit),
              names_from=`Baseline SARS-CoV-2`, 
              values_from=c(N, Responder, `GMT/GMC`)) %>% 
  inner_join(rrdiff_bl, c("Arm", "subgroup", "Marker", "Visit")) %>% 
  inner_join(rgmt_bl, c("Arm", "Marker", "Visit")) %>% 
  select(Arm, Visit, Marker, `N_Positive`, `Responder_Positive`, 
         `GMT/GMC_Positive`, `N_Negative`, `Responder_Negative`, 
         `GMT/GMC_Negative`, rrdiff, `Ratios of GMT/GMC`) %>% 
  arrange(Arm, Visit, Marker)

tab_vacc <- tab_bl %>% dplyr::filter(Arm == "Vaccine") %>% select(-Arm)
tab_plcb <- tab_bl %>% dplyr::filter(Arm == "Placebo") %>% select(-Arm)

print("Done with table15") 

save(tlf, tab_dm_pos, tab_dm_neg, tab_bind, tab_pseudo, tab_wt, tab_gm,
     tab_gmr, tab_rgmt, tab_rrdiff, tab_neg, tab_pos, tab_vacc, tab_plcb,
     file = here::here("output", "Tables.Rdata"))


