#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(dplyr)
library(tidyr)

# DB: Scheduled for deletion
# library(COVIDcorr)
# data(dat.mock)
dat.mock <- read.csv(here("..", "data_raw", data_name))

# select subset: two phase samples
plot_dat <- subset(dat.mock, TwophasesampInd==1)
saveRDS(plot_dat, file = here("data_clean", "plot_data.rds"))

# wide to long format by marker (bindSpike, bindRBD, ID50, ID80) and time (Baseline, Day29, Day57, Delta29overB, Delta57overB)
# add label to variables: Bserostatus, Trt, Time
# define strata: age >= 65, risk, sex at birth(1=female, 0=male), RaceEthnic, Dich_RaceEthnic, age x risk
plot_dat_long <- plot_dat %>%
  select(Ptid, Trt, Bserostatus, EventIndPrimaryD29, EventIndPrimaryD57, Perprotocol, 
         Age, HighRiskInd, Sex, EthnicityHispanic,
         EthnicityNotreported, EthnicityUnknown, WhiteNonHispanic,
         BbindSpike, Day29bindSpike, Day57bindSpike, Delta29overBbindSpike, Delta57overBbindSpike,
         BbindRBD, Day29bindRBD, Day57bindRBD, Delta29overBbindRBD, Delta57overBbindRBD,
         Bpseudoneutid50, Day29pseudoneutid50, Day57pseudoneutid50, Delta29overBpseudoneutid50, Delta57overBpseudoneutid50,
         Bpseudoneutid80, Day29pseudoneutid80, Day57pseudoneutid80, Delta29overBpseudoneutid80, Delta57overBpseudoneutid80, 
  ) %>%
  pivot_longer(!Ptid:WhiteNonHispanic, names_to = "time_marker", values_to = "value") %>%
  mutate(time = factor(gsub("bindSpike|bindRBD|pseudoneutid50|pseudoneutid80", "", time_marker), levels=c("B", "Day29","Day57", "Delta29overB", "Delta57overB"), labels=c("Day 1", "Day 29","Day 57", "Delta29overB", "Delta57overB")),
         marker = factor(gsub("^B|^Day29|^Day57|^Delta29overB|^Delta57overB", "", time_marker), levels=c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80"), labels=c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")),
         Bserostatus = factor(Bserostatus, levels = c(0, 1), labels = c("Baseline Neg","Baseline Pos")),
         Trt = factor(Trt, levels = c(0, 1), labels = c("Placebo","Vaccine"))
  ) %>%
  mutate(AgeInd = ifelse(Age>=65, "Age>=65","Age<65"),
         HighRiskInd = factor(HighRiskInd, levels = c(0, 1), labels = c("Not at risk","At risk")),
         Sex = factor(Sex, levels = c(0, 1), labels = c("Male","Female")),
         RaceEthnic = ifelse(WhiteNonHispanic==1, "White Non-Hispanic", ifelse(WhiteNonHispanic==0, "Comm. of Color", NA)),
         Dich_RaceEthnic = ifelse(EthnicityHispanic==1, "Hispanic or Latino", ifelse(EthnicityHispanic==0 & EthnicityNotreported==0 & EthnicityUnknown==0, "Not Hispanic or Latino", NA))) %>%
  mutate(AgeInd_HighRiskInd = paste(AgeInd, HighRiskInd))
saveRDS(plot_dat_long, file = here("data_clean", "plot_dat_long.rds"))

# stack Intercurrent cases and Perprotocol cohort
# define event (cases, non-cases):
#     Intercurrent cases: EventIndPrimaryD29==1 & EventIndPrimaryD57==0
#     Per-protocol cases: EventIndPrimaryD29==1 & EventIndPrimaryD57==1
# lower bound cutoff:      LLoQ   0.5*LLoQ
#     bindSpike, bindRBD   34     17
#     pseudoneutid50  	   49     25
#     pseudoneutid80  	   43     22
# upper bound cutoff:      ULoQ
#     bindSpike, bindRBD   19136250
# positive threshold:  
#     bindSpike, bindRBD   34 (LLoQ) 
#     pseudoneutid50, 80   20 (LLoD)
# response: 
#     positive if (baseline < positive threshold and post-baseline > positive threshold) 
#              or (baseline > positive threshold and 4 fold increase at post-baseline) 
#     negative
plot_dat_long_stacked <-  plot_dat_long %>% 
  filter(EventIndPrimaryD29==1 & EventIndPrimaryD57==0) %>% 
  mutate(CohortInd="Intercurrent",
         Event = 1) %>%
  bind_rows(plot_dat_long %>% 
              filter(Perprotocol==1) %>% 
              filter(EventIndPrimaryD29==EventIndPrimaryD57) %>% 
              mutate(CohortInd="PP") %>%
              mutate(Event = ifelse(EventIndPrimaryD29==1 & EventIndPrimaryD57==1, 1, 0))) %>%
  mutate(Event = factor(Event, levels = c(1, 0), labels = c("Cases","Non-cases"))) %>%
  mutate(cohort_event = paste(CohortInd, Event)) %>%
  mutate(HalfLLoQ = ifelse(marker=="pseudoneutid50", log10(25), 
                           ifelse(marker=="pseudoneutid80", log10(22), 
                                  ifelse(marker %in% c("bindSpike","bindRBD"), log10(17), NA))),
         LLoQ = ifelse(marker=="pseudoneutid50", log10(49), 
                       ifelse(marker=="pseudoneutid80", log10(43), 
                              ifelse(marker %in% c("bindSpike","bindRBD"), log10(34), NA))),
         ULoQ = ifelse(marker %in% c("bindSpike","bindRBD"), log10(19136250), NA),
         
         pos_threshold = ifelse(marker %in% c("bindSpike","bindRBD"), log10(34),
                                ifelse(marker %in% c("pseudoneutid50","pseudoneutid80"), log10(20), NA)),
         
         baseline_lt_thres = ifelse(time=="Day 1" & value >= pos_threshold, 1, 0),
         increase_4F_D29 = ifelse(time=="Delta29overB" & value>log10(4), 1, 0), 
         increase_4F_D57 = ifelse(time=="Delta57overB" & value>log10(4), 1, 0)) %>%
  group_by(Ptid, marker) %>%
  mutate(baseline_lt_thres_ptid=max(baseline_lt_thres),
         increase_4F_D29_ptid=max(increase_4F_D29),
         increase_4F_D57_ptid=max(increase_4F_D57)) %>%
  ungroup() %>%
  filter(time %in% c("Day 1","Day 29","Day 57")) %>%
  mutate(response = ifelse(baseline_lt_thres_ptid == 0 & value >= pos_threshold, 1,
                           ifelse(baseline_lt_thres_ptid == 1 & time == "Day 1", 1, 
                                  ifelse(baseline_lt_thres_ptid == 1 & time == "Day 29" & increase_4F_D29_ptid==1, 1, 
                                         ifelse(baseline_lt_thres_ptid == 1 & time == "Day 57" & increase_4F_D57_ptid==1, 1,0)))),
         value2 = ifelse(value<LLoQ, HalfLLoQ, 
                         ifelse(!is.na(ULoQ) & value>ULoQ, ULoQ, value)))  # value2 is the truncated version of value
saveRDS(plot_dat_long_stacked, file = here("data_clean", "plot_dat_long_stacked.rds"))  

#### figure specific data prep
# 1. define response rate:
# 2. make subsample datasets such that the jitter plot for each subgroup in each panel <= 25 data points

#### for Figure 1. intercurrent vs pp, case vs non-case, (Day 1), Day 29 Day 57
groupby_vars1=c("Trt", "Bserostatus", "Event", "CohortInd", "time", "marker")

plot_dat_long_stacked_plot1 <- 
  plot_dat_long_stacked %>% group_by_at(groupby_vars1) %>%
  mutate(num = sum(response), 
         denom=n(), 
         RespRate = paste0(num,"/",denom,"=\n",round(num/denom*100, 1),"%"))
saveRDS(plot_dat_long_stacked_plot1, file = here("data_clean", "plot_dat_long_stacked_plot1.rds"))  


plot.25sample1 <- plot_dat_long_stacked_plot1 %>% 
  group_by_at(groupby_vars1) %>%
  sample_n((ifelse(n()>=25, 25, n())), replace=F) %>% filter(time=="Day 57") %>% 
  ungroup() %>%
  select(c("Ptid", groupby_vars1[!groupby_vars1 %in% "time"])) %>%
  inner_join(plot_dat_long_stacked_plot1, by=c("Ptid", groupby_vars1[!groupby_vars1 %in% "time"]))
saveRDS(plot.25sample1, file = here("data_clean", "plot.25sample1.rds"))  

#### for Figure 3. intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57, by if Age >=65 and if at risk
groupby_vars3 <- c("Trt", "Bserostatus", "Event", "CohortInd", "time", "marker", "AgeInd", "HighRiskInd")

plot_dat_long_stacked_plot3 <- 
  plot_dat_long_stacked %>% group_by_at(groupby_vars3) %>%
  mutate(num = sum(response), 
         denom=n(), 
         RespRate = paste0(num,"/",denom,"=\n",round(num/denom*100, 1),"%"))
saveRDS(plot_dat_long_stacked_plot3, file = here("data_clean", "plot_dat_long_stacked_plot3.rds"))  

plot.25sample3 <-  plot_dat_long_stacked_plot3 %>% 
  group_by_at(groupby_vars3) %>%
  sample_n((ifelse(n()>=25, 25, n())), replace=F) %>% filter(time=="Day 57") %>% 
  ungroup() %>%
  select(c("Ptid", groupby_vars3[!groupby_vars3 %in% "time"])) %>%
  inner_join(plot_dat_long_stacked_plot3, by=c("Ptid", groupby_vars3[!groupby_vars3 %in% "time"]))
saveRDS(plot.25sample3, file = here("data_clean", "plot.25sample3.rds"))