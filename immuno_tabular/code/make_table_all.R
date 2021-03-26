##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here("..")) #
source(here::here("..", "_common.R")) #
##################################################

# Immunogenicity Tables

# Reload clean_data
base::load(here::here("data_clean", "params.Rdata"))
base::load(here::here("data_clean", "ds_all.Rdata"))
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

sub_grp_col <- c("subgroup", "Arm", "Baseline SARS-CoV-2", "Group", 
                 "Ind", "Marker", "Visit")

resp_v <- c(grep("Resp|2llod|4llod", unique(ds_resp_l$resp_cat), value = T) %>%
              grep("bind", ., value = T),
            grep("Resp|FR2|FR4", unique(ds_resp_l$resp_cat), value = T) %>%
              grep("neut", ., value = T))

rpcnt <- ds_resp_l %>% 
  dplyr::filter(Visit != "Day 1" & resp_cat %in% resp_v) %>% 
  group_split(across(all_of(sub_grp_col))) %>%
  map_dfr(get_rr, stratum="tps.stratum", weights="wt.subcohort", sub_grp_col=sub_grp_col)

tab_rr <- inner_join(
  rpcnt,
  ds_resp_l %>%
    mutate(rspndr = response*wt.subcohort) %>% 
    group_by(across(all_of(sub_grp_col))) %>%
    summarise(N = n(), Nw = sum(wt.subcohort), rspndr = sum(rspndr),
              .groups = 'drop'),
  by = sub_grp_col) %>%
  mutate(rslt = case_when(
    is.na(ci_l)|is.na(ci_u) ~ 
      sprintf("%s/%s = %.1f%%", round(rspndr,1), round(Nw,1), response*100),
    TRUE ~ 
      sprintf("%s/%s = %.1f%%\n(%.1f%%, %.1f%%)", 
              round(rspndr, 1), round(Nw, 1), response*100, ci_l*100, ci_u*100))
  ) %>%
  pivot_wider(
    id_cols = c(subgroup, Group, Arm, `Baseline SARS-CoV-2`, Marker, Visit, N),
    names_from = Ind, values_from = rslt) %>% 
  arrange(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker)


tab_bind <- tab_rr %>% 
  dplyr::filter(Marker %in% labels_all$Marker[grep("bind", labels_all$marker)]) %>% 
  select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, N, Responder, 
         `% Greater than 2xLLOD`,`% Greater than 4xLLOD`)

print("Done with table2") 

# Table 3 & 4. For the ID50 pseudo-virus & MN50 WT live virus neutralization 
# antibody marker, the estimated percentage of participants defined as 
# responders, participants with % 2-Fold Rise (2FR), and participants with 4-fold 
# rise (4FR) will be provided with the corresponding 95% CIs 
# 
# Output: tab_pseudo & tab_wt

tab_pseudo <- tab_rr %>% 
  dplyr::filter(Marker  %in% labels_all$Marker[grep("pseudoneutid50", labels_all$marker)]) %>% 
  select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, N, Responder,
         `% 2-Fold Rise`, `% 4-Fold Rise`)
tab_wt <- tab_rr %>% 
  dplyr::filter(Marker  %in% labels_all$Marker[grep("liveneutmn50", labels_all$marker)]) %>% 
  select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, N, Responder,
         `% 2-Fold Rise`, `% 4-Fold Rise`)

print("Done with table3 & 4") 

# Table 5. Geometric mean titers (GMTs) and geometric mean concentrations (GMCs)
# will be summarized along with their 95% CIs using the t-distribution
# approximation of log-transformed concentrations/titers (for each of the 5
# Spike-targeted marker types including pseudovirus-nAb ID50 and ID80
# and WT live virus-nAb MN50, as well as for binding Ab to N).
# 
# Output: tab_gm
gm_v <- c(assays_col, 
          gmr_v <- grep("Delta", unique(ds_mag_l$mag_cat), value = T) %>% 
            grep("overB", ., value = T))

sub_grp_col <- c("subgroup", "Arm", "Baseline SARS-CoV-2", "Group", "mag_cat")

rgm <- ds_mag_l %>% 
  dplyr::filter(mag_cat %in% gm_v) %>% 
  group_split(across(all_of(sub_grp_col))) %>%
  map_dfr(get_gm,weights="wt.subcohort",stratum="tps.stratum",sub_grp_col=sub_grp_col)

tab_gm <- inner_join(
  rgm,
  ds_mag_l %>%
    group_by(subgroup, mag_cat, Group, Visit, Marker, Arm, `Baseline SARS-CoV-2`) %>%
    summarise(N = n(), .groups = 'drop'),
  by = sub_grp_col
) 

tab_gmt <- tab_gm %>% 
  dplyr::filter(!grepl("Delta", mag_cat)) %>% 
  mutate(`GMT/GMC`= sprintf("%.0f\n(%.0f, %.0f)", 10^mag, 10^ci_l, 10^ci_u)) %>% 
  arrange(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker) %>% 
  select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, N, `GMT/GMC`) 

print("Done with table5") 

### Table 6. GMTRs/GMCRs will be summarized with 95% CI (t-distribution 
# approximation) for any post-baseline values compared to baseline, and
# post-Day 57 values compared to Day 57
# 
# Output: tab_gmr

tab_gmt_gmr <- inner_join(
  tab_gmt %>% 
    filter(Visit == "Day 1") %>% 
    select(-Visit) %>% 
    rename(`Baseline GMT/GMC` = `GMT/GMC`),
  tab_gmt %>% 
    filter(Visit != "Day 1") %>% 
    rename(`Post Baseline GMT/GMC` = `GMT/GMC`),
  by = c("subgroup", "Arm", "Group", "Baseline SARS-CoV-2", "N", "Marker")
) %>% 
  mutate(Visit = paste0(gsub("ay ", "", Visit), " fold-rise over D1"))

tab_gmr <- tab_gm %>% 
  dplyr::filter(grepl("Delta", mag_cat)) %>%
  mutate(`GMTR/GMCR`=sprintf("%.2f\n(%.2f, %.2f)", 10^mag, 10^ci_l, 10^ci_u)) %>% 
  select(subgroup, Group, Arm, `Baseline SARS-CoV-2`, Visit, N, Marker, `GMTR/GMCR`) %>% 
  inner_join(tab_gmt_gmr,
             by = c("subgroup", "Group", "Arm", "Baseline SARS-CoV-2", "Visit", "N", "Marker")
  ) %>% 
  arrange(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker) %>% 
  select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, N,  
         `Baseline GMT/GMC`, `Post Baseline GMT/GMC`, `GMTR/GMCR`) %>% 
  rename(`Baseline\nGMT/GMC`=`Baseline GMT/GMC`, 
         `Post Baseline\nGMT/GMC`=`Post Baseline GMT/GMC`,)

print("Done with table6") 

### Table 7. The ratios of GMTs/GMCs will be estimated between groups with the
# two-sided 95% CIs calculated using t-distribution approximation of 
# log-transformed titers/concentrations
# Output: tab_rgmt
# 
# Ratios of GMT/GMC between subgroups among vacinees

ds_mag_l_rgmt <- ds_mag_l %>%
  dplyr::filter(subgroup %in% c("Age", "Risk for Severe Covid-19", 
                                "Age, Risk for Severe Covid-19",  
                                "Sex", 
                                "Hispanic or Latino ethnicity",
                                "Underrepresented minority status") & 
                !grepl("Delta", mag_cat) &
                Group != "Not reported and unknown ") %>%
  mutate(
    subgroup = case_when(
      Group %in% c("Age $<$ 65 At-risk", "Age $<$ 65 Not at-risk") ~ "Age, Risk 1",
      Group %in% c("Age $\\geq$ 65 At-risk", "Age $\\geq$ 65 Not at-risk") ~ "Age, Risk 2",
      TRUE ~ as.character(subgroup)),
    comp_i = as.character(subgroup), 
    Group = factor(Group, 
                   levels = c("Age $\\geq$ 65", "Age $<$ 65", 
                              "At-risk", "Not at-risk",    
                              "Age $<$ 65 At-risk", "Age $<$ 65 Not at-risk", 
                              "Age $\\geq$ 65 At-risk", "Age $\\geq$ 65 Not at-risk", 
                              "Male", "Female", 
                              "Hispanic or Latino", "Not Hispanic or Latino",  
                              "Communities of Color", "White Non-Hispanic")))

comp_v <- "Group"
f_v <- as.formula(sprintf("mag ~ %s", comp_v))
sub_grp_col <- c("subgroup", "Arm", "Baseline SARS-CoV-2", "Visit", "Marker")

rgmt <- ds_mag_l_rgmt %>%
  group_split(across(all_of(sub_grp_col))) %>%
  map_dfr(get_rgmt, comp_v=comp_v, f_v = f_v, sub_grp_col = sub_grp_col, 
          comp_lev=NULL, stratum = "tps.stratum", weights = "wt.subcohort") 

tab_gmtx <- tab_gmt %>% 
  dplyr::filter(subgroup %in% c("Age", 
                                "Risk for Severe Covid-19", 
                                "Age, Risk for Severe Covid-19",  
                                "Sex", 
                                "Hispanic or Latino ethnicity",
                                "Underrepresented minority status") &
                  Group != "Not reported and unknown ") %>%
  mutate(subgroup = case_when(
    Group %in% c("Age $<$ 65 At-risk", "Age $<$ 65 Not at-risk") ~ 
      "Age, Risk 1",
    Group %in% c("Age $\\geq$ 65 At-risk", "Age $\\geq$ 65 Not at-risk") ~ 
      "Age, Risk 2",
    TRUE ~ as.character(subgroup)),
    Group = factor(Group, 
                   levels = c("Age $\\geq$ 65", "Age $<$ 65", 
                              "At-risk", "Not at-risk",    
                              "Age $<$ 65 At-risk", "Age $<$ 65 Not at-risk", 
                              "Age $\\geq$ 65 At-risk", "Age $\\geq$ 65 Not at-risk", 
                              "Male", "Female", 
                              "Hispanic or Latino", "Not Hispanic or Latino",  
                              "Communities of Color", "White Non-Hispanic")))
  
tab_gmt_rgmt <- tab_gmtx %>% 
  mutate(groupn = 2-match(Group, unique(sort(Group)))%%2) %>% 
  pivot_wider(id_cols = c(subgroup, Arm, `Baseline SARS-CoV-2`, Visit, Marker),
              names_from = groupn, values_from = `GMT/GMC`, 
              names_prefix = "Group")

tab_rgmt <- inner_join(rgmt, 
             tab_gmt_rgmt, 
             by = c("Baseline SARS-CoV-2", "Arm", "subgroup", "Marker", "Visit")) %>% 
  rename(`Group 1 vs 2` = comp, 
         `Group 1 GMT/GMC` = `Group1`, 
         `Group 2 GMT/GMC` = `Group2`) %>% 
  mutate(`Ratios of GMT/GMC` = sprintf("%.2f\n(%.2f, %.2f)", 
                                       10^Estimate, 10^ci_l, 10^ci_u),
         subgroup = factor(subgroup, levels = c("Age",  
                                                "Risk for Severe Covid-19",
                                                "Age, Risk 1", "Age, Risk 2",
                                                "Sex",
                                                "Hispanic or Latino ethnicity",
                                                "Underrepresented minority status"))) %>% 
  select(`Group 1 vs 2`, subgroup, Visit, Arm, `Baseline SARS-CoV-2`, Marker, 
         `Group 1 GMT/GMC`, `Group 2 GMT/GMC`, `Ratios of GMT/GMC`) %>% 
  arrange(subgroup, Visit, Arm, `Baseline SARS-CoV-2`, Marker)

print("Done with table7") 

### Table 8. The differences in the responder rates, 2FRs, 4FRs between groups 
# will be computed along with the two-sided 95% CIs by the Wilson-Score
# method without continuity correction (Newcombe, 1998).
# Output: tab_rrdiff

comp_v <- c("Vaccine", "Placebo")
sub_grp_col <- c("subgroup", "Arm", "Baseline SARS-CoV-2", "Group", 
                  "Ind", "Marker", "Visit")
rpcnt_diff <- ds_resp_l %>% 
  dplyr::filter(Visit != "Day 1" & Ind %in% c("Responder", "% 2-Fold Rise", "% 4-Fold Rise") &
                subgroup=="All participants") %>% 
  group_split(across(all_of(sub_grp_col))) %>%
  map_dfr(get_rr, weights="wt.subcohort", stratum="tps.stratum", sub_grp_col=sub_grp_col)

rrdiff <- rpcnt_diff %>% 
  dplyr::filter(subgroup == "All participants") %>% 
  mutate(Arm = match(as.character(Arm), comp_v), comp = paste(comp_v, collapse = " vs ")) %>% 
  pivot_wider(names_from = Arm, values_from = c(response, ci_l, ci_u), names_sep = "") %>% 
  mutate(Estimate = response1-response2,
         ci_l = Estimate-sqrt((response1-ci_l1)^2+(response2-ci_u2)^2),
         ci_u = Estimate+sqrt((response1-ci_u1)^2+(response2-ci_l2)^2)) %>% 
  select(-c(response1, response2, ci_l1, ci_l2, ci_u1, ci_u2))

tab_rrdiff <- rrdiff %>% 
  mutate(rslt = sprintf("%s%%\n(%s%%, %s%%)", round(Estimate * 100, 1), 
                        round(ci_l*100, 1), round(ci_u*100, 1)),
    Comparison = comp
  ) %>% 
  pivot_wider(-c(Estimate, ci_u, ci_l),
              names_from = Ind,
              values_from = rslt) %>%
  arrange(subgroup, Group, Visit, `Baseline SARS-CoV-2`, Marker, Comparison) %>% 
  select(subgroup, Group, Visit, `Baseline SARS-CoV-2`, Marker, Comparison, Responder,
         `% 2-Fold Rise`, `% 4-Fold Rise`)

print("Done with table8") 

# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021)

# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021) 
# Vaccine vs Placebo within Baseline status
ds_mag_l_Rx <- ds_mag_l %>%
  dplyr::filter(subgroup == "All participants" & !grepl("Delta", mag_cat))

comp_v <- "Arm"
# contrasts(ds_mag_l_Rx$Arm) <- contr.treatment(2, base = 2)
f_v <- as.formula(sprintf("mag ~ %s", comp_v))
sub_grp_col_Rx <- c("subgroup", "Group", "Baseline SARS-CoV-2", "mag_cat", "Marker", "Visit")

rgmt_Rx <- ds_mag_l_Rx %>%
  group_split(across(all_of(sub_grp_col_Rx))) %>%
  map_dfr(get_rgmt, comp_v=comp_v, comp_lev=c("Vaccine", "Placebo"),
          f_v = f_v, sub_grp_col = sub_grp_col_Rx, 
          stratum = "tps.stratum", weights = "wt.subcohort") 

tab_rgmt_Rx <- rgmt_Rx %>%
  mutate(`Ratios of GMT/GMC` = sprintf("%.2f\n(%.2f, %.2f)", 
                                       10^Estimate, 10^ci_l, 10^ci_u)) %>% 
  select(`Baseline SARS-CoV-2`, Marker, Visit, subgroup, Group, comp, `Ratios of GMT/GMC`) %>% 
  arrange(subgroup, `Baseline SARS-CoV-2`, Visit, Marker)

tab_key_Rx <- full_join(tab_rr, 
                        tab_gmt,
                        by = c("subgroup", "Arm", "Baseline SARS-CoV-2", 
                               "Group", "Visit", "N", "Marker")) %>%
  dplyr::filter(as.character(subgroup) == "All participants") %>%
  pivot_wider(id_cols = c(subgroup, Group, `Baseline SARS-CoV-2`, Marker, Visit),  
              names_from = Arm, 
              values_from = c(N, Responder, `GMT/GMC`)
  ) %>% 
  # Join with RR difference (tab_rr)
  inner_join(tab_rrdiff, 
             by = c("subgroup", "Group", "Baseline SARS-CoV-2", "Visit", "Marker")
  ) %>% 
  # Join with GMT/GMC Ratios (tab_rgmt) 
  inner_join(tab_rgmt_Rx, 
             by = c("subgroup", "Group", "Baseline SARS-CoV-2", "Visit", "Marker")
  ) %>% 
  select(`Baseline SARS-CoV-2`, Visit, Marker, `N_Vaccine`, `Responder_Vaccine`, 
         `GMT/GMC_Vaccine`, `N_Placebo`, `Responder_Placebo`, `GMT/GMC_Placebo`, 
         Responder, `Ratios of GMT/GMC`) %>% 
  arrange(`Baseline SARS-CoV-2`, Visit, Marker)

tab_neg <- filter(tab_key_Rx, `Baseline SARS-CoV-2` == "Negative") %>% 
  select(-c(`Baseline SARS-CoV-2`))
tab_pos <- filter(tab_key_Rx, `Baseline SARS-CoV-2` == "Positive") %>% 
  select(-c(`Baseline SARS-CoV-2`))

print("Done with table13-14") 

###################################################

comp_v <- c("Positive", "Negative")

rrdiff_bl <- rpcnt_diff %>% 
  dplyr::filter(subgroup == "All participants" & Ind == "Responder") %>% 
  mutate(`Baseline SARS-CoV-2` = match(as.character(`Baseline SARS-CoV-2`), comp_v), 
         comp = paste(comp_v, collapse = " vs ")) %>% 
  pivot_wider(names_from = `Baseline SARS-CoV-2`, 
              values_from = c(response, ci_l, ci_u), names_sep = "") %>% 
  mutate(Estimate = response1-response2,
         ci_l = Estimate-sqrt((response1-ci_l1)^2+(response2-ci_u2)^2),
         ci_u = Estimate+sqrt((response1-ci_u1)^2+(response2-ci_l2)^2)) %>% 
  select(-c(response1, response2, ci_l1, ci_l2, ci_u1, ci_u2))

tab_rrdiff_bl <- rrdiff_bl %>% 
  mutate(rslt = case_when(
    !complete.cases(ci_l, ci_u) ~ sprintf("%s%%", round(Estimate * 100, 1)),
    complete.cases(ci_l, ci_u) ~ sprintf("%s%%\n(%s%%, %s%%)", round(Estimate * 100, 1), 
                   round(ci_l*100, 1), round(ci_u*100, 1)))
  ) %>% 
  select(-c(Estimate, ci_u, ci_l)) 


ds_mag_l_bl <- ds_mag_l %>%
  dplyr::filter(subgroup == "All participants" & !grepl("Delta", mag_cat))

contrasts(ds_mag_l_bl$`Baseline SARS-CoV-2`) <- contr.treatment(2, base = 2)
# f_v <- as.formula("mag ~ `Baseline SARS-CoV-2`")
sub_grp_col_bl <- c("subgroup", "Arm", "Visit", "Marker", "mag_cat")

rgmt_bl <- ds_mag_l_bl %>%
  group_split(across(all_of(sub_grp_col_bl))) %>%
  map_dfr(get_rgmt, 
          comp_v = "Baseline SARS-CoV-2", comp_lev = c("Positive", "Negative"),
          f_v = "mag ~ `Baseline SARS-CoV-2`", sub_grp_col = sub_grp_col_bl, 
          stratum = "tps.stratum", weights = "wt.subcohort") 

tab_rgmt_bl <- rgmt_bl %>%
  mutate(`Ratios of GMT/GMC` = sprintf("%.2f\n(%.2f, %.2f)", 
                                       10^Estimate, 10^ci_l, 10^ci_u)) %>% 
  select(Arm, Marker, Visit, subgroup, comp, `Ratios of GMT/GMC`) %>% 
  arrange(subgroup, Arm, Visit, Marker)

 tab_bl <- full_join(tab_rr, 
                     tab_gmt, 
                     by = c("Arm", "subgroup", "Group",
                             "Baseline SARS-CoV-2", "Visit", "Marker", "N")) %>%
  dplyr::filter(as.character(subgroup) == "All participants") %>%
  pivot_wider(id_cols = c(Arm, subgroup, Group, `Baseline SARS-CoV-2`, Marker, Visit),
              names_from=`Baseline SARS-CoV-2`, 
              values_from=c(N, Responder, `GMT/GMC`)
  ) %>% 
  # Join with RR difference (tab_rr)
  inner_join(tab_rrdiff_bl, 
             c("Arm", "subgroup", "Group", "Marker", "Visit")
  ) %>% 
  # Join with GMT/GMC Ratios (tab_rgmt) 
  inner_join(tab_rgmt_bl, 
             by = c("Arm", "subgroup", "Marker", "Visit", "comp")
  ) %>% 
  select(Arm, Visit, Marker, `N_Positive`, `Responder_Positive`, 
         `GMT/GMC_Positive`, `N_Negative`, `Responder_Negative`, 
         `GMT/GMC_Negative`, rslt, `Ratios of GMT/GMC`) %>% 
  arrange(Arm, Visit, Marker)

tab_vacc <- tab_bl %>% dplyr::filter(Arm == "Vaccine") %>% select(-Arm)
tab_plcb <- tab_bl %>% dplyr::filter(Arm == "Placebo") %>% select(-Arm)

print("Done with table15") 

save(tlf, tab_dm_pos, tab_dm_neg, tab_bind, tab_pseudo, tab_wt, tab_gmt,
     tab_gmr, tab_rgmt, tab_rrdiff, tab_neg, tab_pos, tab_vacc, tab_plcb,
     file = here::here("output", "Tables.Rdata"))


