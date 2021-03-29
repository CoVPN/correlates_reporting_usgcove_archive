###########################
# File: cor_graphical_verification_code.R
# Author: Di Lu
# Creation Date: 2021/03/22
# Last Edit Date: 2021/03/28
# Description: Verification of the plots targeted for
# independent double programming
# NEWS:
#
#
#
#

#renv::activate()

# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(PResiduals)
library(mosaic)
library(psych)
library(data.table)
########################################################################################################

# initial data processing
#dat <- read.csv("D:/Verification_Di/correlates_reporting/immuno_graphical/verification/practice_data.csv")

dat <- read.csv(
  here("cor_graphical/verification/verification_input/practice_data.csv")
)



dat <- dat %>% 
  mutate(cohort_event = ifelse(EventIndPrimaryD29 == 1 & EventIndPrimaryD57 == 0, "Intercurrent Cases",
                               ifelse(Perprotocol == 1 & EventIndPrimaryD29 == 1 & EventIndPrimaryD57 == 1, "PP Cases",
                                      ifelse(Perprotocol == 1 & EventIndPrimaryD29 == 0 & EventIndPrimaryD57 == 0, "PP Non-cases", NA)
  )
))


dat <- dat %>% drop_na(cohort_event)


# wide to long
data.variables_needed_bindSpike <- dat %>%
  select(
    Ptid, 
    Trt, 
    MinorityInd, 
    EthnicityHispanic, 
    EthnicityNotreported, 
    EthnicityUnknown, 
    HighRiskInd, 
    Age, 
    BMI, 
    Sex,
    Bserostatus, 
    Fullvaccine, 
    Perprotocol, 
    EventIndPrimaryD29, 
    EventIndPrimaryD57,
    SubcohortInd, 
    age.geq.65, 
    TwophasesampInd,
    Bstratum, 
    wt, 
    wt.2, 
    race, 
    WhiteNonHispanic, 
    cohort_event,
    BbindSpike, 
    Day29bindSpike, 
    Day57bindSpike, 
    Delta29overBbindSpike, 
    Delta57overBbindSpike, 
    Delta57over29bindSpike
  ) %>%
  mutate(assay = "bindSpike") %>%
  rename(
    B = BbindSpike, 
    Day29 = Day29bindSpike,
    Day57 = Day57bindSpike, 
    Delta29overB = Delta29overBbindSpike, 
    Delta57overB = Delta57overBbindSpike, 
    Delta57over29 = Delta57over29bindSpike)


data.variables_needed_bindRBD <- dat %>%
  select(
    Ptid, 
    Trt, 
    MinorityInd, 
    EthnicityHispanic, 
    EthnicityNotreported, 
    EthnicityUnknown, 
    HighRiskInd, 
    Age, 
    BMI, 
    Sex,
    Bserostatus, 
    Fullvaccine, 
    Perprotocol, 
    EventIndPrimaryD29, 
    EventIndPrimaryD57, 
    SubcohortInd, 
    age.geq.65, 
    TwophasesampInd,
    Bstratum, 
    wt, 
    wt.2, 
    race, 
    WhiteNonHispanic, 
    cohort_event,
    BbindRBD, 
    Day29bindRBD, 
    Day57bindRBD, 
    Delta29overBbindRBD, 
    Delta57overBbindRBD, 
    Delta57over29bindRBD
  ) %>%
  mutate(assay = "bindRBD") %>%
  rename(
    B = BbindRBD, 
    Day29 = Day29bindRBD, 
    Day57 = Day57bindRBD, 
    Delta29overB = Delta29overBbindRBD, 
    Delta57overB = Delta57overBbindRBD, 
    Delta57over29 = Delta57over29bindRBD)


data.variables_needed_pseudoneutid50 <- dat %>%
  select(
    Ptid, 
    Trt, 
    MinorityInd, 
    EthnicityHispanic, 
    EthnicityNotreported, 
    EthnicityUnknown, 
    HighRiskInd, 
    Age, 
    BMI, 
    Sex,
    Bserostatus, 
    Fullvaccine, 
    Perprotocol, 
    EventIndPrimaryD29, 
    EventIndPrimaryD57, 
    SubcohortInd, 
    age.geq.65, 
    TwophasesampInd,
    Bstratum, 
    wt, 
    wt.2, 
    race, 
    WhiteNonHispanic, 
    cohort_event,
    Bpseudoneutid50, 
    Day29pseudoneutid50, 
    Day57pseudoneutid50, 
    Delta29overBpseudoneutid50, 
    Delta57overBpseudoneutid50, 
    Delta57over29pseudoneutid50
  ) %>%
  mutate(assay = "pseudoneutid50") %>%
  rename(
    B = Bpseudoneutid50, 
    Day29 = Day29pseudoneutid50, 
    Day57 = Day57pseudoneutid50, 
    Delta29overB = Delta29overBpseudoneutid50, 
    Delta57overB = Delta57overBpseudoneutid50, 
    Delta57over29 = Delta57over29pseudoneutid50)


data.variables_needed_pseudoneutid80 <- dat %>%
  select(
    Ptid, 
    Trt, 
    MinorityInd, 
    EthnicityHispanic, 
    EthnicityNotreported, 
    EthnicityUnknown, 
    HighRiskInd, 
    Age, 
    BMI, 
    Sex,
    Bserostatus, 
    Fullvaccine, 
    Perprotocol, 
    EventIndPrimaryD29, 
    EventIndPrimaryD57, 
    SubcohortInd, 
    age.geq.65, 
    TwophasesampInd,
    Bstratum, 
    wt, 
    wt.2, 
    race, 
    WhiteNonHispanic, 
    cohort_event,
    Bpseudoneutid80, 
    Day29pseudoneutid80, 
    Day57pseudoneutid80, 
    Delta29overBpseudoneutid80, 
    Delta57overBpseudoneutid80, 
    Delta57over29pseudoneutid80
  ) %>%
  mutate(assay = "pseudoneutid80") %>%
  rename(
    B = Bpseudoneutid80, 
    Day29 = Day29pseudoneutid80, 
    Day57 = Day57pseudoneutid80, 
    Delta29overB = Delta29overBpseudoneutid80, 
    Delta57overB = Delta57overBpseudoneutid80, 
    Delta57over29 = Delta57over29pseudoneutid80)

dat.long <- rbind(
  data.variables_needed_bindSpike, 
  data.variables_needed_bindRBD, 
  data.variables_needed_pseudoneutid50, 
  data.variables_needed_pseudoneutid80)

dat.cor.subset <- dat %>% 
  filter(TwophasesampInd == 1)

dat.long.cor.subset <- dat.long %>% 
  filter(TwophasesampInd == 1)

dat.long.cor.subset <- dat.long.cor.subset %>% 
  mutate(Dich_RaceEthnic = ifelse(EthnicityHispanic == 1, "Hispanic or Latino",
  ifelse(EthnicityHispanic == 0 & EthnicityNotreported == 0 & EthnicityUnknown == 0, "Not Hispanic or Latino", NA)
))


dat.long.cor.subset <- dat.long.cor.subset %>% mutate(LLoD = ifelse(assay == "pseudoneutid50", log10(10),
  ifelse(assay == "pseudoneutid80", log10(10),
    ifelse(assay == "bindSpike", log10(20),
      ifelse(assay == "bindRBD", log10(20), NA)
    )
  )
))



dat.long.cor.subset <- dat.long.cor.subset %>%
  group_by(assay) %>%
  summarise(ULoQ = max(Day29, Day29)) %>%
  inner_join(dat.long.cor.subset)


dat.long.cor.subset <- dat.long.cor.subset %>% mutate(Delta29overB = Day29 - B, Delta57overB = Day57 - B)

dat.long.cor.subset <- dat.long.cor.subset %>%
  mutate(demo_lab = paste(dat.long.cor.subset$age.geq.65, dat.long.cor.subset$HighRiskInd, sep = " ")) %>%
  mutate(demo_lab = as.factor(demo_lab))

dat.long.cor.subset <- dat.long.cor.subset %>%
  mutate(trt_bstatus_label = paste(dat.long.cor.subset$Trt, dat.long.cor.subset$Bserostatus, sep = " ")) %>%
  mutate(trt_bstatus_label = as.factor(trt_bstatus_label))

dat.long.cor.subset <- dat.long.cor.subset %>%
  mutate(age_geq_65_label = ifelse(age.geq.65 >= 65, ">=65", "<65")) %>%
  mutate(highrisk_label = ifelse(HighRiskInd == 1, "at risk", "not at risk")) %>%
  mutate(sex_label = ifelse(Sex == 1, "Male", "Female")) %>%
  mutate(
    ethnicity_label = ifelse(EthnicityHispanic == 1, "Hispanic or Latino", 
                             ifelse(EthnicityHispanic == 0 & EthnicityNotreported == 0 & EthnicityUnknown == 0, "Not Hispanic or Latino", "Not reported and unknown"))) %>%
  mutate(minority_label = ifelse(WhiteNonHispanic == 1, "White Non-Hispanic", "Comm. of Color"))

dat.long.cor.subset <- dat.long.cor.subset %>%
  mutate(age_risk_label = paste(dat.long.cor.subset$age_geq_65_label, dat.long.cor.subset$highrisk_label, sep = " ")) %>%
  mutate(age_sex_label = paste(dat.long.cor.subset$age_geq_65_label, dat.long.cor.subset$sex_label, sep = " ")) %>%
  mutate(age_minority_label = paste(dat.long.cor.subset$age_geq_65_label, dat.long.cor.subset$minority_label, sep = " "))


dat.longer.cor.subset <- dat.long.cor.subset %>% select(
  Ptid, Trt, Bserostatus, EventIndPrimaryD29, EventIndPrimaryD57, Perprotocol,
  cohort_event, Age, age_geq_65_label, highrisk_label, age_risk_label, sex_label,
  minority_label, Dich_RaceEthnic, assay, LLoD, wt, wt.2, B, Day29, Day57, Delta29overB, Delta57overB
)


dat.longer.cor.subset <- gather(dat.longer.cor.subset, time, value, B:Delta57overB, factor_key = TRUE) %>%
  mutate(time = case_when(
    time == "B" ~ "Day 1",
    time == "Day29" ~ "Day 29",
    time == "Day57" ~ "Day 57",
    time == "Delta29overB" ~ "Delta29overB",
    time == "Delta57overB" ~ "Delta57overB"
  ))


dat.longer.cor.subset <- dat.longer.cor.subset %>%
  mutate(baseline_lt_thres = ifelse(time == "Day 1" & value >= LLoD, 1, 0)) %>%
  mutate(increase_4F_D29 = ifelse(time == "Delta29overB" & value > log10(4), 1, 0)) %>%
  mutate(increase_4F_D57 = ifelse(time == "Delta57overB" & value > log10(4), 1, 0))


dat.longer.cor.subset <- dat.longer.cor.subset %>%
  group_by(Ptid, assay) %>%
  mutate(
    baseline_lt_thres_ptid = max(baseline_lt_thres, na.rm = T),
    increase_4F_D29_ptid = max(increase_4F_D29, na.rm = T),
    increase_4F_D57_ptid = max(increase_4F_D57, na.rm = T)
  ) %>%
  ungroup()


dat.longer.cor.subset <- dat.longer.cor.subset %>% filter(time == "Day 1" | time == "Day 29" | time == "Day 57")

dat.longer.cor.subset <- dat.longer.cor.subset %>% 
  mutate(
  response = ifelse((baseline_lt_thres_ptid == 0 & value >= LLoD) |
  (baseline_lt_thres_ptid == 1 & time == "Day 1") |
  (baseline_lt_thres_ptid == 1 & time == "Day 29" & increase_4F_D29_ptid == 1) |
  (baseline_lt_thres_ptid == 1 & time == "Day 57" & increase_4F_D57_ptid == 1), 1, 0))




#lineplots of Binding Antibody to Spike: baseline negative vaccine arm (2 timepoints)


dat.longer.cor.subset.plot1 <- readRDS(
  here("cor_graphical/verification/verification_input/longer_cor_data_plot1.rds")
)

plot.25sample1 <- readRDS(
  here("cor_graphical/verification/verification_input/longer_cor_data_plot1.rds")
)


lineplots_neg_vaccine_bindSpike <- dat.longer.cor.subset.plot1 %>% 
  filter(Trt == "Vaccine", 
         Bserostatus == "Baseline Neg", 
         assay == "bindSpike")

lineplots_neg_vaccine_bindSpike <- lineplots_neg_vaccine_bindSpike %>% 
  filter(time == "Day 29" | time == "Day 57")

p <- ggplot(lineplots_neg_vaccine_bindSpike, aes(x = time, y = value)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), limits = c(0, 10)) +
  facet_wrap(~cohort_event)
#geom_line(aes(x=time, y=value,group=Ptid))

lineplots_neg_vaccine_bindSpike %>% 
  select(cohort_event,
         time,
         RespRate) %>% 
  distinct()

#violinplots of Pseudovirus Neutralization ID50: baseline negative vaccine arm (2 timepoints)

lineplots_neg_vaccine_pseudoneutid50 <- dat.longer.cor.subset.plot1 %>% 
  filter(Trt == "Vaccine", 
         Bserostatus == "Baseline Neg", 
         assay == "pseudoneutid50")

lineplots_neg_vaccine_pseudoneutid50 <- lineplots_neg_vaccine_pseudoneutid50 %>% 
  filter(time == "Day 29" | time == "Day 57")

p <- ggplot(lineplots_neg_vaccine_pseudoneutid50, aes(x = time, y = value)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), limits = c(0, 10)) +
  facet_wrap(~cohort_event)
# geom_line(aes(x=time, y=value,group=Ptid))
p

lineplots_neg_vaccine_pseudoneutid50 %>% 
  select(cohort_event,
         time,
         RespRate) %>% 
  distinct()


#lineplots of Binding Antibody to Spike: baseline negative vaccine arm (3 timepoints)
lineplots_neg_vaccine_bindSpike <- dat.longer.cor.subset.plot1 %>% filter(Trt == "Vaccine", Bserostatus == "Baseline Neg", assay == "bindSpike")

lineplots_neg_vaccine_bindSpike <- lineplots_neg_vaccine_bindSpike

p <- ggplot(lineplots_neg_vaccine_bindSpike, aes(x = time, y = value)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), limits = c(0, 10)) +
  facet_wrap(~cohort_event)
# geom_line(aes(x=time, y=value,group=Ptid))
p

lineplots_neg_vaccine_bindSpike %>% 
  select(cohort_event,
         time,
         RespRate) %>% 
  distinct()





