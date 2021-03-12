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
library(COVIDcorr)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

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

# Stack a Rx = "Total" to the original data
ds_long_ttl <- bind_rows(
  ds_long %>% mutate(Rx = "Total"),
  ds_long
)

# Calculate % for categorical covariates
dm_cat <- inner_join(
  ds_long_ttl %>%
    group_by(Rx, subgroup, subgroup_cat) %>%
    summarise(n = n(), .groups = 'drop'),
  ds_long_ttl %>%
    group_by(Rx, subgroup) %>%
    summarise(N = n(), .groups = 'drop'),
  by = c("Rx", "subgroup")
) %>%
  mutate(pct = n / N,
         rslt1 = sprintf("%s (%.1f%%)", n, n / N * 100), 
         rslt2 = sprintf("%s/%s = %.1f%%", n, N, n / N * 100)) %>% 
  dplyr::filter(subgroup %in% cat_v)

# Calculate mean and range for numeric covariates
dm_num <- ds_long_ttl %>%
  distinct_at(all_of(c("Rx", "Ptid", num_v1, num_v2))) %>%
  pivot_longer(cols = all_of(c(num_v1, num_v2)), 
               names_to = "subgroup", 
               values_to = "subgroup_cat") %>%
  group_by(Rx, subgroup) %>%
  summarise(
    min = min(subgroup_cat, na.rm = T), 
    max = max(subgroup_cat, na.rm = T),
    mean = mean(subgroup_cat, na.rm = T),
    sd = sd(subgroup_cat, na.rm = T), 
    rslt1 = sprintf("%.1f (%.1f, %.1f)", mean, min, max),
    rslt2 = sprintf("%.1f +/- %.1f", mean, sd),
    N = n(),
    .groups = 'drop'
  ) %>% 
  mutate(subgroup_cat = case_when(subgroup %in% num_v1 ~ "Mean (Range)",
                                  subgroup %in% num_v2 ~ "Mean +/- SD"))

tab_dm <- full_join(
  dm_cat, 
  dm_num,
  by = c("Rx", "subgroup", "subgroup_cat","N", "rslt1", "rslt2")
) %>%
  mutate(rslt = case_when(subgroup %in% cat_v ~ rslt1,
                          subgroup %in% num_v1 ~ rslt1,
                          subgroup %in% num_v2 ~ rslt2),
         subgroup = factor(subgroup, 
                           levels = c("Age", "BMI", "Sex", "Race", 
                                      "Hispanic or Latino ethnicity", 
                                      "Risk for Severe Covid-19")),
         Rx = factor(Rx, levels = c("Placebo", "Vaccine", "Total"))) %>% 
  pivot_wider(c(Rx, subgroup, subgroup_cat, rslt),
              names_from = Rx, 
              names_sort = T,
              values_from = rslt) %>%
  arrange(subgroup, subgroup_cat)

names(tab_dm)<- c("subgroup", "Characteristics",
                  sprintf("Placebo\n(N = %s)", sum(ds_s$Trt==0)),
                  sprintf("Vaccine\n(N = %s)", sum(ds_s$Trt==1)),
                  sprintf("Total\n(N = %s)", nrow(ds_s)))

### Table 2. Responder Rates & Proportions of Magnitudes >= 2FR, 4FR
# For each binding antibody marker, the estimated percentage of participants
# defined as responders, and with concentrations >= 2x LLOD or >=
# 4 x LLOD, will be provided with the corresponding 95% CIs
# 
# Output: tab_bind


# Variables used for stratification in the tables
# subgroup: SAP Table 6: Baseline Subgroups
# Rx and Baseline: Assigned treatment Arms * Baseline COVID-19 Status
# Group: Category in each subgroup

sub_grp_col <- c("subgroup", "Rx", "Baseline", "Group", "resp_cat")
sub_grp <- as.formula(paste0("~", paste(sub_grp_col, collapse = "+")))

resp_v <- c(grep("Resp|2llod|4llod", unique(ds_resp_l$resp_cat), value = T) %>% 
              grep("bind", ., value = T),
            grep("Resp|FR2|FR4", unique(ds_resp_l$resp_cat), value = T) %>% 
              grep("neut", ., value = T))

dssvy <- svydesign(ids = ~ Ptid, 
                   # strata = ~ Wstratum, 
                   weights = ~ wt.subcohort,
                   data = ds_resp_l %>% filter(resp_cat %in% resp_v &
                                                 Visit != "Day 1"),
                   nest = T)

rpcnt <- svyby(~ response, by=sub_grp, dssvy, svyciprop, vartype = "ci")

tab_rr <- inner_join(
  rpcnt,
  ds_resp_l %>%
    mutate(rspndr = response*wt.subcohort) %>% 
    group_by(Rx, Baseline, Marker, resp_cat, subgroup, Group, Visit, Ind) %>%
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
    id_cols = c(subgroup, Group, Rx, Baseline, Marker, Visit, N),
    names_from = Ind, values_from = rslt) 

tab_bind <- tab_rr %>% 
  dplyr::filter(Marker %in% c("Anti N IgG (IU/ml)", "Anti RBD IgG (IU/ml)", 
                              "Anti Spike IgG (IU/ml)")) %>% 
  select(Rx, Baseline, subgroup, Group, Marker, Visit, N, Responder, 
         `% Greater than 2xLLOD`,`% Greater than 4xLLOD`)

# Table 3 & 4. For the ID50 pseudo-virus & MN50 WT live virus neutralization 
# antibody marker, the estimated percentage of participants defined as 
# responders, participants with % 2-Fold Rise (2FR), and participants with 4-fold 
# rise (4FR) will be provided with the corresponding 95% CIs 
# 
# Output: tab_pseudo & tab_wt

tab_pseudo <- tab_rr %>% dplyr::filter(Marker == "Pseudovirus-nAb ID50") %>% 
  select(Rx, Baseline, subgroup, Group, Marker, Visit, N, Responder,
         `% 2-Fold Rise`, `% 4-Fold Rise`)
tab_wt <- tab_rr %>% dplyr::filter(Marker == "Live virus-nAb MN50") %>% 
  select(Rx, Baseline, subgroup, Group, Marker, Visit, N, Responder,
         `% 2-Fold Rise`, `% 4-Fold Rise`)


# Table 5. Geometric mean titers (GMTs) and geometric mean concentrations (GMCs)
# will be summarized along with their 95% CIs using the t-distribution
# approximation of log-transformed concentrations/titers (for each of the 5
# Spike-targeted marker types including pseudovirus-nAb ID50 and ID80
# and WT live virus-nAb MN50, as well as for binding Ab to N).
# 
# Output: tab_gm
gm_v <- c(bAb_v, pnAb_v, lnAb_v, 
          gmr_v <- grep("DeltaDay", names(ds), value = T, fixed = F) %>% 
            grep("overB", ., value = T))

sub_grp_col <- c("subgroup", "Rx", "Baseline", "Group", "mag_cat")
sub_grp <- as.formula(paste0("~", paste(sub_grp_col, collapse = "+")))

dssvy <- svydesign(ids = ~ Ptid, 
                   # strata = ~ Wstratum, 
                   weights = ~ wt.subcohort,
                   data = ds_mag_l %>% filter(mag_cat %in% gm_v),
                   nest = T)

rgm <- svyby(~ mag, sub_grp, dssvy, svymean, vartype = "ci")

tab_gm <- inner_join(
  rgm,
  ds_mag_l %>%
    group_by(subgroup, mag_cat, Group, Visit, Marker, Rx, Baseline) %>%
    summarise(N = n(), .groups = 'drop'),
  by = sub_grp_col
) 

tab_gmt <- tab_gm %>% 
  dplyr::filter(!grepl("Delta", mag_cat)) %>% 
  mutate(`GMT/GMC`= sprintf("%.0f\n(%.0f, %.0f)", 10^mag, 10^ci_l, 10^ci_u)) %>% 
  select(subgroup, Rx, Group, Baseline, Visit, N, Marker, `GMT/GMC`) 

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
  by = c("subgroup", "Rx", "Group", "Baseline", "N", "Marker")
) %>% 
  mutate(Visit = paste0(gsub("ay ", "", Visit), " fold-rise over D1"))


tab_gmr <- tab_gm %>% 
  dplyr::filter(grepl("Delta", mag_cat)) %>%
  mutate(`GMTR/GMCR`=sprintf("%.2f\n(%.2f, %.2f)", 10^mag, 10^ci_l, 10^ci_u)) %>% 
  select(subgroup, Group, Rx, Baseline, Visit, N, Marker, `GMTR/GMCR`) %>% 
  inner_join(tab_gmt_gmr,
             by = c("subgroup", "Group", "Rx", "Baseline", "Visit", "N", "Marker")
  ) %>% 
  select(subgroup, Group, Rx, Baseline, Visit, N, Marker, 
         `Baseline GMT/GMC`, `Post Baseline GMT/GMC`, `GMTR/GMCR`)

### Table 7. The ratios of GMTs/GMCs will be estimated between groups with the
# two-sided 95% CIs calculated using t-distribution approximation of 
# log-transformed titers/concentrations
# Output: tab_gmtrC
# 
# Ratios of GMT/GMC between subgroups among vacinees
ds_mag_l_gmtr <- ds_mag_l %>%
  dplyr::filter(subgroup %in% c("Age",  "Sex", 
                                "Risk for Severe Covid-19", 
                                "Race and ethnic group") & 
                  !grepl("Delta", mag_cat)) %>% 
  mutate(comp_i = as.character(subgroup), 
         subgroup = factor(comp_i, levels = c("Age", "Sex", 
                                              "Risk for Severe Covid-19",
                                              "Hispanic or Latino ethnicity",
                                              "Race and ethnic group")),
         Group = factor(Group, 
                        levels = c("Age >= 65", "Age < 65", "At-risk", 
                                   "Not at-risk", "Female", "Male",
                                   "Communities of Color",
                                   "White Non-Hispanic")))

comp_v <- "Group"
f_v <- as.formula(sprintf("mag ~ %s", comp_v))
sub_grp_col <- c("subgroup", "Rx", "Baseline", "Visit", "Marker")

gmtr <- ds_mag_l_gmtr %>%
  group_split(across(all_of(sub_grp_col))) %>%
  map_dfr(get_gmtr, comp_v=comp_v, f_v = f_v, sub_grp_col = sub_grp_col, 
          weights = "wt.subcohort") 


tab_gmt_gmtr <- inner_join(
  tab_gmt, 
  tab_gmt %>% 
    group_by(subgroup) %>% 
    summarise(Group = unique(Group),
              ngroup = n_distinct(Group), 
              groupn = match(Group, unique(sort(Group)))),
  by = c("subgroup", "Group")
) %>%
  dplyr::filter(ngroup == 2) %>% 
  pivot_wider(id_cols = c(subgroup, Rx, Baseline, Visit, Marker),
              names_from = groupn, values_from = `GMT/GMC`, 
              names_prefix = "Group")

tab_gmtr <- gmtr %>%
  mutate(ci_u = 10^(Estimate+1.96*`Std. Error`), 
         ci_l = 10^(Estimate-1.96*`Std. Error`),
         `Ratios of GMT/GMC` = sprintf("%.2f\n(%.2f, %.2f)", 
                                       10^Estimate, ci_l, ci_u)) %>% 
  inner_join(tab_gmt_gmtr, 
             by = c("Baseline", "Rx", "subgroup", "Marker", "Visit")) %>% 
  rename(Comparison = comp, 
         `Group 1 GMT/GMC` = `Group1`, 
         `Group 2 GMT/GMC` = `Group2`) %>% 
  select(subgroup, Rx, Baseline, Marker, Visit, Comparison, `Group 1 GMT/GMC`,
         `Group 2 GMT/GMC`, `Ratios of GMT/GMC`) %>% 
  arrange(subgroup, Rx, Baseline, Visit, Marker)

### Table 8. The differences in the responder rates, 2FRs, 4FRs between groups 
# will be computed along with the two-sided 95% CIs by the Wilson-Score
# method without continuity correction (Newcombe, 1998).
# Output: tab_rrdiff
# 
comp_v <- "Rx"
f_v <- as.formula(sprintf("response ~ %s", comp_v))
contrasts(ds_resp_l$Rx) <- contr.treatment(2, base = 2)
sub_grp_col_diff <- c("Baseline", "subgroup", "Group", "Ind", "Visit", "Marker")

rrdiff <- ds_resp_l %>%
  dplyr::filter(Visit != "Day 1") %>% 
  group_split(across(all_of(sub_grp_col_diff))) %>%
  map_dfr(get_rrdiff, comp_v = comp_v, f_v = f_v, weights = "wt.subcohort", 
          sub_grp_col = sub_grp_col_diff)

tab_rrdiff <- rrdiff %>% 
  mutate(
    ci_u = Estimate + 1.96 * `Std. Error`, 
    ci_l = Estimate - 1.96 * `Std. Error`,
    rslt = sprintf("%s%%\n(%s%%, %s%%)", round(Estimate * 100, 1), 
                   round(ci_l*100, 1), round(ci_u*100, 1))
  ) %>% 
  pivot_wider(-c(Estimate, `Std. Error`, ci_u, ci_l),
              names_from = Ind,
              values_from = rslt) %>%
  select(subgroup, Group, Baseline, Visit, Marker, comp, Responder,
         `% 2-Fold Rise`, `% 4-Fold Rise`)

# 8b Responder rate differences between cases vs non-cases & 95% CI of 
# Titers or Concentrations
comp_v <- "Group"
f_v <- as.formula(sprintf("response ~ %s", comp_v))
sub_grp_col_case <- c("Rx", "Baseline","subgroup", "Ind", "Marker", "Visit")

ds_resp_l_case <- ds_resp_l %>% 
  dplyr::filter(subgroup == "Event" & Ind == "Responder") %>% 
  mutate(Group = factor(Group, levels = c("Cases", "Non-Cases")))

contrasts(ds_resp_l_case$Group) <- contr.treatment(2, base = 2)

rrdiff_case <- ds_resp_l_case %>%
  group_split(across(all_of(sub_grp_col_case))) %>%
  map_dfr(get_rrdiff, comp_v = comp_v, f_v = f_v, weights = "weight", 
          sub_grp_col = sub_grp_col_case)

tab_rrdiff_case <- rrdiff_case %>% 
  mutate(
    ci_u = Estimate + 1.96 * `Std. Error`, 
    ci_l = Estimate - 1.96 * `Std. Error`,
    rslt = sprintf("%s%%\n(%s%%, %s%%)", round(Estimate * 100, 1), 
                   round(ci_l*100, 1), round(ci_u*100, 1))
  ) %>% 
  pivot_wider(-c(Estimate, `Std. Error`, ci_u, ci_l),
              names_from = Ind,
              values_from = rslt) 


# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021)
# Cases vs Non-cases

ds_resp_case <- ds_resp_l %>% 
  dplyr::filter(subgroup == "Event" & grepl("Resp", resp_cat))

resp_v <- unique(ds_resp_case$resp_cat)

dssvy_case <- svydesign(ids = ~ Ptid, 
                        # strata = ~ Wstratum, 
                        weights = ~ weight,
                        data = ds_resp_case,
                        nest = T)

sub_grp_col <- c("subgroup", "Rx", "Baseline", "Group", "resp_cat")
sub_grp <- as.formula(paste0("~", paste(sub_grp_col, collapse = "+")))

rpcnt_case <- svyby(~response, by=sub_grp, dssvy_case, svyciprop, vartype = "ci")

tab_rr_case <- rpcnt_case %>%
  inner_join(
    ds_resp_case %>%
      mutate(rspndr = response*wt.subcohort) %>% 
      group_by(Rx, Baseline, Marker, resp_cat, subgroup, Group, Visit) %>%
      summarise(N = n(), Nw = sum(wt.subcohort), rspndr = sum(rspndr),
                .groups = 'drop'),
    by = c("Rx", "Baseline", "subgroup", "Group","resp_cat")
  ) %>%
  mutate(Responder = case_when(
    is.na(ci_l)|is.na(ci_u) ~ 
      sprintf("%s/%s = %.1f%%", round(rspndr,1), round(Nw,1), response*100),
    TRUE ~ 
      sprintf("%s/%s = %.1f%%\n(%.1f%%, %.1f%%)", 
              round(rspndr, 1), round(Nw, 1), response*100, ci_l*100, ci_u*100))
  ) %>% 
  select(subgroup, Rx, Baseline, Group, Visit, N, Marker, Responder)

#########

gm_v <- c(bAb_v, pnAb_v, lnAb_v)
ds_mag_case <- ds_mag_l %>% 
  dplyr::filter(subgroup == "Event" & mag_cat %in% gm_v & !is.na(weight))

sub_grp_col <- c("subgroup", "Rx", "Baseline", "Group", "mag_cat")
sub_grp <- as.formula(paste0("~", paste(sub_grp_col, collapse = "+")))

dssvy_case <- svydesign(ids = ~ Ptid, 
                        # strata = ~ Wstratum, 
                        weights = ~ weight,
                        data = ds_mag_case,
                        nest = T)

gm_f <- paste0("~", paste(gm_v, collapse = " + "))
rgm_case <- svyby(~ mag, sub_grp, dssvy_case, svymean, vartype = "ci")

tab_gm_case <- rgm_case %>%
  inner_join(ds_mag_case %>%
               group_by(subgroup, mag_cat, Group, Visit, Marker, Rx, Baseline) %>%
               summarise(N = n(), .groups = 'drop'),
             by = c("subgroup", "Group", "Rx", "Baseline", "mag_cat")
  ) %>%
  mutate(`GMT/GMC` = sprintf("%.0f\n(%.0f, %.0f)", 10^mag, 10^ci_l, 10^ci_u)
  ) %>% 
  select(subgroup, Rx, Baseline, Group, Visit, N, Marker, `GMT/GMC`)

###

###
ds_mag_l_gmtr_case <- ds_mag_l %>%
  dplyr::filter(subgroup == "Event" & !is.na(weight))

comp_v <- "Group"
f_v <- as.formula(sprintf("mag ~ %s", comp_v))
sub_grp_col <- c("subgroup", "Rx", "Baseline", "Visit", "Marker")

gmtr_case <- ds_mag_l_gmtr_case %>%
  group_split(across(all_of(sub_grp_col))) %>%
  map_dfr(get_gmtr, comp_v=comp_v, f_v = f_v, sub_grp_col = sub_grp_col, 
          weights = "weight") 

tab_gmtr_case <- gmtr_case %>%
  mutate(ci_u = 10^(Estimate+1.96*`Std. Error`), 
         ci_l = 10^(Estimate-1.96*`Std. Error`),
         `Ratios of GMT/GMC` = sprintf("%.2f\n(%.2f, %.2f)", 
                                       10^Estimate, ci_l, ci_u)) %>% 
  select(Rx, Baseline, Marker, Visit, subgroup, comp, `Ratios of GMT/GMC`) %>% 
  arrange(subgroup, Rx, Baseline, Visit, Marker)

tab_key_case <- full_join(tab_rr_case, 
                          tab_gm_case,
                          by = c("subgroup", "Rx", "Baseline", 
                                 "Group", "N", "Marker", "Visit")) %>%
  pivot_wider(names_from = Group, 
              values_from = c(N, Responder, `GMT/GMC`)) %>% 
  # Join with RR difference (tab_rr)
  full_join(tab_rrdiff_case, 
            by = c("Rx", "subgroup", "Baseline", "Visit", "Marker")
  ) %>% 
  # Join with GMT/GMC Ratios (tab_gmtr) 
  full_join(tab_gmtr_case, 
            by = c("Rx", "Baseline", "Visit", "Marker", "subgroup")
  ) %>% 
  select(Rx, Baseline, Visit, Marker, `N_Non-Cases`, `Responder_Non-Cases`, 
         `GMT/GMC_Non-Cases`, `N_Cases`, `Responder_Cases`, `GMT/GMC_Cases`, 
         Responder, `Ratios of GMT/GMC`) %>% 
  arrange(Rx, Baseline) 

case_vacc_neg <- filter(tab_key_case, Rx == "Vaccine", Baseline == "Negative") %>% 
  select(-c(Rx, Baseline))
case_vacc_pos <- filter(tab_key_case, Rx == "Vaccine", Baseline == "Positive") %>% 
  select(-c(Rx, Baseline))
case_plcb_pos <- filter(tab_key_case, Rx == "Placebo", Baseline == "Positive") %>% 
  select(-c(Rx, Baseline))


# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021) 
# Vaccine vs Placebo within Baseline status
ds_mag_l_Rx <- ds_mag_l %>%
  dplyr::filter(subgroup == "All participants" & !grepl("Delta", mag_cat))

comp_v <- "Rx"
contrasts(ds_mag_l_Rx$Rx) <- contr.treatment(2, base = 2)
f_v <- as.formula(sprintf("mag ~ %s", comp_v))
sub_grp_col_Rx <- c("subgroup", "Group", "Baseline", "mag_cat", "Marker", "Visit")

gmtr_Rx <- ds_mag_l_Rx %>%
  group_split(across(all_of(sub_grp_col_Rx))) %>%
  map_dfr(get_gmtr, comp_v=comp_v, f_v = f_v, sub_grp_col = sub_grp_col_Rx, 
          weights = "wt.subcohort", desc = F) 

tab_gmtr_Rx <- gmtr_Rx %>%
  mutate(ci_u = 10^(Estimate+1.96*`Std. Error`), 
         ci_l = 10^(Estimate-1.96*`Std. Error`),
         `Ratios of GMT/GMC` = sprintf("%.2f\n(%.2f, %.2f)", 
                                       10^Estimate, ci_l, ci_u)) %>% 
  select(Baseline, Marker, Visit, subgroup, Group, comp, `Ratios of GMT/GMC`) %>% 
  arrange(subgroup, Baseline, Visit, Marker)

tab_key_Rx <- full_join(tab_rr, 
                        tab_gmt,
                        by = c("subgroup", "Rx", "Baseline", 
                               "Group", "Visit", "N", "Marker")) %>%
  dplyr::filter(as.character(subgroup) == "All participants") %>%
  pivot_wider(id_cols = c(subgroup, Group, Baseline, Marker, Visit),  
              names_from = Rx, 
              values_from = c(N, Responder, `GMT/GMC`)
  ) %>% 
  # Join with RR difference (tab_rr)
  inner_join(tab_rrdiff, 
             by = c("subgroup", "Group", "Baseline", "Visit", "Marker")
  ) %>% 
  # Join with GMT/GMC Ratios (tab_gmtrC) 
  inner_join(tab_gmtr_Rx, 
             by = c("subgroup", "Group", "Baseline", "Visit", "Marker")
  ) %>% 
  select(Baseline, Visit, Marker, `N_Vaccine`, `Responder_Vaccine`, 
         `GMT/GMC_Vaccine`, `N_Placebo`, `Responder_Placebo`, `GMT/GMC_Placebo`, 
         Responder, `Ratios of GMT/GMC`) %>% 
  arrange(Baseline, Visit, Marker)

tab_neg <- filter(tab_key_Rx, Baseline == "Negative") %>% 
  select(-c(Baseline))
tab_pos <- filter(tab_key_Rx, Baseline == "Positive") %>% 
  select(-c(Baseline))

###################################################

comp_v <- "Baseline"
f_v <- as.formula(sprintf("response ~ %s", comp_v))
sub_grp_col_bl <- c("Rx", "subgroup", "Marker", "Visit", "resp_cat")

ds_resp_l_bl <- ds_resp_l %>% 
  dplyr::filter(subgroup == "All participants" & Ind == "Responder")

rrdiff_bl <- ds_resp_l_bl %>%
  group_split(across(all_of(sub_grp_col_bl))) %>%
  map_dfr(get_rrdiff, comp_v = comp_v, f_v = f_v, weights = "wt.subcohort", 
          sub_grp_col = sub_grp_col_bl) 

tab_rrdiff_bl <- rrdiff_bl %>% 
  mutate(
    ci_u = Estimate + 1.96 * `Std. Error`, 
    ci_l = Estimate - 1.96 * `Std. Error`,
    rslt = sprintf("%s%%\n(%s%%, %s%%)", round(Estimate * 100, 1), 
                   round(ci_l*100, 1), round(ci_u*100, 1))
  ) %>% 
  select(-c(Estimate, `Std. Error`, ci_u, ci_l)) 


ds_mag_l_bl <- ds_mag_l %>%
  dplyr::filter(subgroup == "All participants" & 
                  Rx == "Vaccine" & !grepl("Delta", mag_cat))

comp_v <- "Baseline"
contrasts(ds_mag_l_bl$Baseline) <- contr.treatment(2, base = 2)
f_v <- as.formula(sprintf("mag ~ %s", comp_v))
sub_grp_col_bl <- c("subgroup", "Group", "Rx", "Visit", "Marker")


gmtr_bl <- ds_mag_l_bl %>%
  group_split(across(all_of(sub_grp_col_bl))) %>%
  map_dfr(get_gmtr, comp_v=comp_v, f_v = f_v, sub_grp_col = sub_grp_col_bl, 
          weights = "wt.subcohort", desc = F) 

tab_gmtr_bl <- gmtr_bl %>%
  mutate( ci_u = 10^(Estimate+1.96*`Std. Error`), 
          ci_l = 10^(Estimate-1.96*`Std. Error`),
          `Ratios of GMT/GMC` = sprintf("%.2f\n(%.2f, %.2f)", 
                                        10^Estimate, ci_l, ci_u)) %>% 
  select(Rx, Marker, Visit, subgroup, comp, `Ratios of GMT/GMC`) %>% 
  arrange(subgroup, Rx, Visit, Marker)

 tab_bl <- full_join(tab_rr, 
                     tab_gmt, 
                     by = c("Rx", "subgroup", "Group",
                             "Baseline", "Visit", "Marker", "N")) %>%
  dplyr::filter(as.character(subgroup) == "All participants") %>%
  pivot_wider(id_cols = c(Rx, subgroup, Group, Baseline, Marker, Visit),
              names_from=Baseline, 
              values_from=c(N, Responder, `GMT/GMC`)
  ) %>% 
  # Join with RR difference (tab_rr)
  inner_join(tab_rrdiff_bl, 
             by = c("subgroup", "Rx","Visit", "Marker")
  ) %>% 
  # Join with GMT/GMC Ratios (tab_gmtrC) 
  inner_join(tab_gmtr_bl, 
             by = c("Visit", "Rx","Marker", "subgroup")
  ) %>% 
  select(Rx, Visit, Marker, `N_Negative`, `Responder_Negative`, 
         `GMT/GMC_Negative`, `N_Positive`, `Responder_Positive`, 
         `GMT/GMC_Positive`, rslt, `Ratios of GMT/GMC`) %>% 
  arrange(Visit, Marker)

tab_vacc <- tab_bl %>% dplyr::filter(Rx == "Vaccine")
tab_plcb <- tab_bl %>% dplyr::filter(Rx == "Placebo")



save(tlf, tab_dm, tab_bind, tab_pseudo, tab_wt, tab_gmt,
     tab_gmr, tab_gmtr, tab_rrdiff, case_vacc_neg, case_vacc_pos, case_plcb_pos, 
     tab_neg, tab_pos, tab_vacc, tab_plcb,
     file = here::here("output", "Tables.Rdata"))

