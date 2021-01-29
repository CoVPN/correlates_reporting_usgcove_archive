##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here(".."))       #
source(here::here("..", "_common.R"))            #
##################################################

# Immunogenicity Tables

# Reload clean_data
base::load(here::here("data_clean", "params.Rdata"))
base::load(here::here("data_clean", "ds_all.Rdata"))

library(survey)
library(tidyverse)

# Select the covariates to be summarised. 
# num_v are columns from ds_long; 
# cat_v are rows of `subgroup`
num_v <- c("Age", "BMI")
cat_v <- c("Age", "Risk for Severe Covid-19", "Sex", "Hispanic or Latino ethnicity", "Race", "Risk for Severe Covid-19")

# Stack a Rx="Total" to the original data
ds_long_ttl <- bind_rows(
  ds_long %>% mutate(Rx="Total"), 
  ds_long
) 

# Calculate % for categorical covariates
dm_cat <- inner_join(
  ds_long_ttl %>% 
    filter(subgroup %in% cat_v) %>% 
    group_by(Rx, subgroup, subgroup_cat) %>% 
    summarise(n=n()),
  
  ds_long_ttl %>% 
    filter(subgroup %in% cat_v) %>% 
    group_by(Rx, subgroup) %>% 
    summarise(N=n())
) %>% 
  mutate(pct=n/N, rslt=sprintf("%s (%.1f)", n, n/N*100), rslt2=sprintf("%s/%s=%.1f%%", n, N, n/N*100))

# Calculate mean and range for numeric covariates
# Mean (Min, Max)
dm_num1 <- ds_long_ttl %>%
  distinct_at(all_of(c("Rx", "Ptid", num_v))) %>% 
  pivot_longer(cols=all_of(num_v), names_to="subgroup", values_to="subgroup_cat") %>% 
  group_by(Rx, subgroup) %>% 
  summarise(min=round(min(subgroup_cat), 1), 
            max=round(max(subgroup_cat), 1),
            mean=round(mean(subgroup_cat), 1), 
            range=sprintf("(%s, %s)", min, max),
            rslt=paste(mean, range),
            N=n(),
            subgroup_cat="Mean (range)")

# Mean +/- SD
dm_num2 <- ds_long_ttl %>%
  distinct_at(all_of(c("Rx", "Ptid", num_v))) %>% 
  pivot_longer(cols=all_of(num_v), names_to="subgroup", values_to="subgroup_cat") %>% 
  group_by(Rx, subgroup) %>% 
  summarise(mean=round(mean(subgroup_cat), 1), 
            sd=round(sd(subgroup_cat), 1),
            rslt=paste(mean, "+/-", sd),
            N=n(),
            subgroup_cat="Mean +/- SD")

tab_dm <- full_join(dm_cat, dm_num1) %>% 
  full_join(dm_num2) %>% 
  pivot_wider(c(Rx, subgroup, subgroup_cat, rslt), names_from=Rx, values_from=rslt) %>% 
  arrange(subgroup, subgroup_cat) %>% 
  relocate(subgroup, subgroup_cat, Placebo, Vaccine, Total) %>% 
  filter(!((subgroup=="Age" & subgroup_cat=="Mean +/- SD")|(subgroup=="BMI" & subgroup_cat=="Mean (range)"))) %>%
  mutate(subgroup = factor(subgroup, levels=c("Age", "BMI", "Sex", "Hispanic or Latino ethnicity", "Race", "Risk for Severe Covid-19"))) %>% 
  arrange(subgroup, subgroup_cat) %>% 
  ungroup()

# Variables used for stratification in the tables
# subgroup: SAP Table 6: Baseline Subgroups
# Group: Category in each subgroup
# Baseline: Baseline COVID-19 Positive or Negative

sub_grp_col <- c("subgroup", "Group", "Baseline")
sub_grp <- as.formula(paste0("~", paste(sub_grp_col, collapse="+"))) 
dssvy <- svydesign(ids=~Bstratum, strata=sub_grp, weights=~wt, data=ds, nest=T)

# 1 - 3 tab_rr: Responder Rates and 95% CI
resp_v <- grep("Resp|FR2|FR4", names(ds), value=T)  
resp_f <- paste0("~", paste(resp_v, collapse = " + "))
rpcnt <- svyby(as.formula(resp_f), sub_grp, dssvy, svymean, vartype="ci")

tab_rr <- rpcnt %>% 
  rename_with(function(x)paste0("rr.",x), !contains(sub_grp_col) & !starts_with("ci_")) %>% 
  pivot_longer(cols=-all_of(sub_grp_col), names_to="resp_cat", values_to = "value") %>% 
  mutate(param=sapply(resp_cat, function(x)strsplit(x, ".", fixed=T)[[1]][1]), 
         resp_cat=sapply(resp_cat, function(x)strsplit(x, ".", fixed=T)[[1]][2])) %>%  
  # Merge with the labels
  inner_join(ds_resp_l %>% 
               group_by(Baseline, Endpoint, resp_cat, subgroup, Group, Visit, Ind) %>%
               summarise(N=n(), rspndr=sum(response)),
             by = c("subgroup", "Group", "Baseline", "resp_cat")) %>%
  # Format into different respond rate columns
  pivot_wider(names_from = param, values_from = value) %>% 
  mutate(rspr= sprintf("%.1f%%", rr*100), ci=sprintf("(%.1f%%, %.1f%%)", ci_l*100, ci_u*100)) %>% 
  pivot_longer(cols=c(rspr, ci), values_to = "rslt") %>% 
  pivot_wider(c(subgroup, Group, Baseline, Visit, Endpoint, name, rslt, N), names_from=Ind, values_from = rslt)

# 4 tab_gm: Geometric Mean & 95% CI of Titers or Concentrations 
gm_f <- paste0("~", paste(c(bAb_v, pnAb_v, lnAb_v), collapse = " + "))
rgm <- svyby(as.formula(gm_f), sub_grp, dssvy, svymean, vartype="ci")

tab_gm <- rgm %>% 
  rename_with(function(x)paste0("gm.",x), !contains(sub_grp_col) & !starts_with("ci_")) %>% 
  pivot_longer(cols=-all_of(sub_grp_col), names_to="mag_cat", values_to = "value") %>% 
  mutate(param=sapply(mag_cat, function(x)strsplit(x, ".", fixed=T)[[1]][1]), 
         mag_cat=sapply(mag_cat, function(x)strsplit(x, ".", fixed=T)[[1]][2])) %>%  
  inner_join(ds_mag_l %>% 
               group_by(subgroup, mag_cat, Group, Visit, Endpoint, Baseline) %>%
               summarise(N=n()),
             by=c("subgroup", "Group", "Baseline", "mag_cat")) %>% 
  pivot_wider(names_from = param, values_from = value) %>% 
  mutate(`GMT/GMC`=sprintf("%.0f", 10^gm), 
         `95% CI`=sprintf("(%.0f, %.0f)", 10^ci_l, 10^ci_u)) 

# 5 - 6: 
gmr_f <- paste0("~", paste(grep("Delta", names(ds), value=T, fixed=F), collapse = " + "))
rgmr <- svyby(as.formula(gmr_f), sub_grp, dssvy, svymean, vartype="ci")


tab_gmr <- rgmr %>% 
  rename_with(function(x)paste0("gmr.",x), !contains(sub_grp_col) & !starts_with("ci_")) %>% 
  pivot_longer(cols=-all_of(sub_grp_col), names_to="colname", values_to = "value") %>% 
  mutate(param = sapply(colname, function(x)strsplit(x, ".", fixed=T)[[1]][1]), 
         colname = sapply(colname, function(x)strsplit(x, ".", fixed=T)[[1]][2])) %>% 
  inner_join(labels_all %>% distinct(colname, Endpoint, Visit, label.long), by="colname") %>% 
  inner_join(ds %>% group_by(across(all_of(sub_grp_col))) %>% summarise(N=n()), by=sub_grp_col) %>% 
  pivot_wider(names_from = param, values_from = value) %>% 
  mutate(`GMTR/GMCR`=sprintf("%.2f", 10^gmr), `95% CI`=sprintf("(%.2f, %.2f)", 10^ci_l, 10^ci_u)) 


# 7a
# Ratios of GMT/GMC between Vaccine vs Placebo 
comp_v <- "Rx"
f_v <- as.formula(sprintf("mag ~ %s", comp_v))

tab_gmtrA <- ds_mag_l %>% 
  filter(subgroup=="All participants" & Visit != "Day 1") %>% 
  group_split(across(c(all_of(sub_grp_col), "Visit", "mag_cat"))) %>% 
  map_dfr(function(x){
    ret <- x %>% 
      group_by(subgroup, Baseline, Endpoint, Visit, mag_cat) %>% 
      summarise(N = n(), 
                estimate = summary(lm(f_v, weights = wt, data=.))$coefficients[2, "Estimate"],
                se=summary(lm(f_v, weights = wt, data=.))$coefficients[2, "Std. Error"],
                ci_u=10^(estimate+1.96*se), ci_l=10^(estimate-1.96*se), 
                `Ratios of GMT/GMC`=round(10^estimate,2),
                `95% CI`=sprintf("(%.2f, %.2f)", ci_l, ci_u))
  })


# 7b Ratios of GMT/GMC between baseline negative vs positive among vacinees
comp_v <- "Baseline"
f_v <- as.formula(sprintf("mag ~ %s", comp_v))
sub_grp_colB <- c("subgroup", "Group", "Rx")

tab_gmtrB <- ds_mag_l %>% 
  filter(subgroup=="All participants" & Rx=="Vaccine" & Visit != "Day 1") %>% 
  group_split(across(c(all_of(sub_grp_colB), "mag_cat"))) %>%
  map_dfr(function(x){
    ret <- x %>% 
      group_by(subgroup, Endpoint, Visit, mag_cat) %>% 
      summarise(N=n(),
                estimate=summary(lm(f_v, weights = wt, data=.))$coefficients[2, "Estimate"],
                se=summary(lm(f_v, weights = wt, data=.))$coefficients[2, "Std. Error"],
                ci_u=10^(estimate+1.96*se), ci_l=10^(estimate-1.96*se), 
                `Ratios of GMT/GMC`=round(10^estimate,2),
                `95% CI`=sprintf("(%.2f, %.2f)", ci_l, ci_u))
  }) 


# 7c
ds_mag_l_7c <- ds_mag_l %>% 
  filter(subgroup %in% c("Age", "Risk for Severe Covid-19", "Age, Risk for Severe Covid-19 ", "Sex", 
                         "Hispanic or Latino ethnicity", "Race and ethnic group")) %>% 
  mutate(comp_i=case_when(as.character(Group) %in% c("Age >= 65 At-risk", "Age >= 65 Not at-risk") ~ "Age >= 65, Risk",
                          as.character(Group) %in% c("Age < 65 At-risk", "Age < 65 Not at-risk") ~ "Age < 65, Risk",
                          TRUE ~ as.character(subgroup)))

comp_v <- "Group"
f_v <- as.formula(sprintf("mag ~ %s", comp_v))
sub_grp_colC <- c("subgroup", "Baseline")

tab_gmtrC <- ds_mag_l_7c %>% 
  filter(Rx=="Vaccine" & Visit != "Day 1") %>% 
  group_split(across(c(all_of(sub_grp_colC), "mag_cat", "comp_i"))) %>% 
  map_dfr(function(x){
    ret <- x %>% 
      group_by(Baseline, Endpoint, mag_cat, Visit, comp_i) %>% 
      summarise(N=n(),
                estimate=summary(lm(f_v, weights = wt, data=.))$coefficients[2, "Estimate"],
                se=summary(lm(f_v, weights = wt, data=.))$coefficients[2, "Std. Error"],
                ci_u=10^(estimate+1.96*se), ci_l=10^(estimate-1.96*se),
                `Ratios of GMT/GMC`=round(10^estimate,2),
                `95% CI`=sprintf("(%.2f, %.2f)", ci_l, ci_u))
  })

tmp <- tab_gmtrC %>% 
  arrange(Baseline, Endpoint, comp_i)
# 8
comp_v <- "Rx"
f_v <- as.formula(sprintf("response ~ %s", comp_v))
contrasts(ds_resp_l$Rx) <- contr.treatment(2, base=2)


tab_rrdiff <- ds_resp_l %>% 
  group_split(across(c(all_of(sub_grp_col), "resp_cat", "Ind"))) %>% 
  map_dfr(function(x){
    ret <- x %>% 
      group_by(Baseline, Endpoint, subgroup, Group, Visit, Ind) %>%
      summarise(N=n(),
                estimate=summary(glm(f_v, weights = wt, data=.))$coefficients[2, "Estimate"],
                se=summary(glm(f_v, weights = wt, data=.))$coefficients[2, "Std. Error"],
                ci_u=estimate+1.96*se, ci_l=estimate-1.96*se,
                est=sprintf("%.1f%%", estimate*100),
                ci=sprintf("(%.1f%%, %.1f%%)", ci_l*100, ci_u*100))}) %>% 
  pivot_longer(cols=c(est, ci), values_to="rslt") %>% 
  pivot_wider(-c(estimate, se, ci_u, ci_l), names_from=Ind, values_from = rslt) 


save(ds_s, tab_dm, tab_rr, tab_gm, tab_gmr, tab_gmtrA, tab_gmtrB, tab_gmtrC, tab_rrdiff, 
     file=here::here("output","Tables.Rdata"))

