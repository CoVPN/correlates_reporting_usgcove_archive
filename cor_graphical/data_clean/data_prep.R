# select subset: two phase samples
plot_dat <- subset(dat.mock, TwophasesampInd==1)

# set flag for intercurrent cases: EventIndPrimaryD29==1 & EventIndPrimaryD57==0
plot_dat$IntercurrentCase <- with(plot_dat,
                                 ifelse(EventIndPrimaryD29==1 & EventIndPrimaryD57==0, 1, 0))

# wide to long format by marker (bindSpike, bindRBD, ID50, ID80) and time (Baseline, Day29, Day57)
# add label to variables: Bserostatus, Trt, Time
# define strata: age >= 65, risk, sex at birth(1=female, 0=male), RaceEthnic, Dich_RaceEthnic, age x risk
plot_dat_long <- plot_dat %>%
  select(Ptid, Trt, Bserostatus, IntercurrentCase, EventIndPrimaryD29, EventIndPrimaryD57, Perprotocol, 
         Age, HighRiskInd, Sex, ethnicity, WhiteNonHispanic,
         BbindSpike, Day29bindSpike, Day57bindSpike,
         BbindRBD, Day29bindRBD, Day57bindRBD,
         Bpseudoneutid50, Day29pseudoneutid50, Day57pseudoneutid50,
         Bpseudoneutid80, Day29pseudoneutid80, Day57pseudoneutid80,
  ) %>%
  pivot_longer(!Ptid:WhiteNonHispanic, names_to = "time_marker", values_to = "value") %>%
  mutate(time = factor(gsub("bindSpike|bindRBD|pseudoneutid50|pseudoneutid80", "", time_marker), levels=c("B", "Day29","Day57"), labels=c("Day 1", "Day 29","Day 57")),
         marker = factor(gsub("^B|^Day29|^Day57", "", time_marker), levels=c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80"), labels=c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")),
         Bserostatus = factor(Bserostatus, levels = c(0, 1), labels = c("Baseline Neg","Baseline Pos")),
         Trt = factor(Trt, levels = c(0, 1), labels = c("Placebo","Vaccine"))
  ) %>%
  mutate(AgeInd = ifelse(Age>=65, "Age>=65","Age<65"),
         HighRiskInd = factor(HighRiskInd, levels = c(0, 1), labels = c("Not at risk","At risk")),
         Sex = factor(Sex, levels = c(0, 1), labels = c("Male","Female")),
         RaceEthnic = factor(ifelse(WhiteNonHispanic==1 & !is.na(WhiteNonHispanic), 1, ifelse(WhiteNonHispanic==0 & !is.na(WhiteNonHispanic), 0, 99))),
         Dich_RaceEthnic = factor(ethnicity)) %>%
  mutate(AgeInd_HighRiskInd = paste(AgeInd, HighRiskInd))


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
plot_dat_long_stacked <-  plot_dat_long %>% 
  filter(IntercurrentCase==1) %>% mutate(CohortInd="Intercurrent") %>% 
  mutate(Event = ifelse(EventIndPrimaryD29==1 & EventIndPrimaryD57==0, 1, 0)) %>%
  bind_rows(plot_dat_long %>% 
              filter(Perprotocol==1) %>% filter(EventIndPrimaryD29==EventIndPrimaryD57) %>% mutate(CohortInd="PP") %>%
              mutate(Event = ifelse(EventIndPrimaryD29==1 & EventIndPrimaryD57==1, 1, 0))) %>%
  mutate(Event = factor(Event, levels = c(1, 0), labels = c("Cases","Non-cases"))) %>%
  mutate(HalfLLoQ = ifelse(marker=="pseudoneutid50", log10(25), 
                           ifelse(marker=="pseudoneutid80", log10(22), 
                                  ifelse(marker %in% c("bindSpike","bindRBD"), log10(17), NA))),
         LLoQ = ifelse(marker=="pseudoneutid50", log10(49), 
                       ifelse(marker=="pseudoneutid80", log10(43), 
                              ifelse(marker %in% c("bindSpike","bindRBD"), log10(34), NA))),
         ULoQ = ifelse(marker %in% c("bindSpike","bindRBD"), log10(19136250), NA)) %>%
  mutate(value2 = ifelse(value<LLoQ, HalfLLoQ, 
                         ifelse(!is.na(ULoQ) & value>ULoQ, ULoQ, value))) %>%  # value2 is the truncated version of value
  mutate(cohort_event = paste(CohortInd, Event))


#### figure specific data prep
# 1. define response rate:
#    binding antibody, a positive response: concentration > 34 IU/ml, a negative response (≤ 34).
#    ID50/80 pseudo/live neut, a positive response: serum ID50 titer > 1:20 (log10(20)), a negative response as the complement
# 2. make subsample datasets such that the jitter plot for each subgroup in each panel <= 25 data points

#### for Figure 1. intercurrent vs pp, case vs non-case, (Day 1), Day 29 Day 57
groupby_vars1=c("Trt", "Bserostatus", "Event", "CohortInd", "time", "marker")

plot_dat_long_stacked_plot1 <- 
  plot_dat_long_stacked %>% group_by_at(groupby_vars1) %>%
  mutate(num = sum(value>=(ifelse(marker %in% c("pseudoneutid50","pseudoneutid80"), log10(20), 
                                  ifelse(marker %in% c("bindSpike","bindRBD"), log10(34), NA)))), 
         denom=n(), 
         RespRate = paste0(num,"/",denom,"=\n",round(num/denom*100, 1),"%"))

plot.25sample1 <- plot_dat_long_stacked_plot1 %>% 
  group_by_at(groupby_vars1) %>%
  sample_n((ifelse(n()>=25, 25, n())), replace=F) %>% filter(time=="Day 57") %>% 
  ungroup() %>%
  select(c("Ptid", groupby_vars1[!groupby_vars1 %in% "time"])) %>%
  inner_join(plot_dat_long_stacked_plot1, by=c("Ptid", groupby_vars1[!groupby_vars1 %in% "time"]))

#### for Figure 3. intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57, by if Age >=65 and if at risk
groupby_vars3 <- c("Trt", "Bserostatus", "Event", "CohortInd", "time", "marker", "AgeInd", "HighRiskInd")

plot_dat_long_stacked_plot3 <- 
  plot_dat_long_stacked %>% group_by_at(groupby_vars3) %>%
  mutate(num = sum(value>=(ifelse(marker %in% c("pseudoneutid50","pseudoneutid80"), log10(20), 
                                  ifelse(marker %in% c("bindSpike","bindRBD"), log10(34), NA)))), 
         denom=n(), 
         RespRate = paste0(num,"/",denom,"=\n",round(num/denom*100, 1),"%"))

plot.25sample3 <-  plot_dat_long_stacked_plot3 %>% 
  group_by_at(groupby_vars3) %>%
  sample_n((ifelse(n()>=25, 25, n())), replace=F) %>% filter(time=="Day 57") %>% 
  ungroup() %>%
  select(c("Ptid", groupby_vars3[!groupby_vars3 %in% "time"])) %>%
  inner_join(plot_dat_long_stacked_plot3, by=c("Ptid", groupby_vars3[!groupby_vars3 %in% "time"]))
