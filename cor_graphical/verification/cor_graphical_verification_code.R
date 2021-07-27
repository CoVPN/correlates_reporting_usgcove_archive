###########################
# File: cor_graphical_verification_code.R
# Author: Di Lu
# Creation Date: 2021/03/22
# Last Edit Date: 2021/07/26
# Description: Verification of the plots targeted for
# independent double programming
# NEWS:
#
#
#
#
renv::activate()

# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(PResiduals)
library(mosaic)
library(psych)
library(data.table)
library(here)
########################################################################################################

# initial data processing


dat <- read.csv(
  here("cor_graphical/verification/data_processed.csv")
)


dat$wt.D57[is.na(dat$wt.D57)] <- 0
dat$wt.D29[is.na(dat$wt.D29)] <- 0


dat <- dat %>% 
  mutate(cohort_event = ifelse(Perprotocol==1 & Bserostatus==0 &  EarlyendpointD29==0 & TwophasesampIndD29==1 & EventIndPrimaryD29==1 & EventTimePrimaryD29 >=7 &  EventTimePrimaryD29 <= (6 + NumberdaysD1toD57 - NumberdaysD1toD29), "Intercurrent Cases",
                               ifelse(Perprotocol==1 & Bserostatus==0 & EarlyendpointD57==0 & TwophasesampIndD57==1 & EventIndPrimaryD57==1, "Post-Peak Cases",
                                      ifelse(Perprotocol==1 &  Bserostatus==0 & EarlyendpointD57==0 & TwophasesampIndD57==1 & EventIndPrimaryD1==0, "Non-Cases", NA)
                               )
  ))


dat <- dat %>% drop_na(cohort_event)


dat1<-dat[,!grepl("*bindN",names(dat))]
dat1<-dat1[,!grepl("*liveneutmn",names(dat1))]
dat1<-dat1[,!grepl("*CPV",names(dat1))]


# wide to long
dat1_bindSpike <- dat1[,!grepl("*bindRBD",names(dat1))]
dat1_bindSpike <- dat1_bindSpike[,!grepl("*pseudoneutid50",names(dat1_bindSpike))]
dat1_bindSpike <- dat1_bindSpike[,!grepl("*pseudoneutid80",names(dat1_bindSpike))]



dat1_bindSpike <- dat1_bindSpike %>%
  mutate(assay = "bindSpike") %>%
  rename(B = BbindSpike,
         Day29 = Day29bindSpike,
         Day57 = Day57bindSpike,
         Delta29overB = Delta29overBbindSpike,
         Delta57overB = Delta57overBbindSpike,
         Delta57over29 = Delta57over29bindSpike
  )


dat1_bindRBD <- dat1[,!grepl("*bindSpike",names(dat1))]
dat1_bindRBD <- dat1_bindRBD[,!grepl("*pseudoneutid50",names(dat1_bindRBD))]
dat1_bindRBD <- dat1_bindRBD[,!grepl("*pseudoneutid80",names(dat1_bindRBD))]



dat1_bindRBD <- dat1_bindRBD %>%
  mutate(assay = "bindRBD") %>%
  rename(B = BbindRBD,
         Day29 = Day29bindRBD,
         Day57 = Day57bindRBD,
         Delta29overB = Delta29overBbindRBD,
         Delta57overB = Delta57overBbindRBD,
         Delta57over29 = Delta57over29bindRBD
  )

dat1_pseudoneutid50 <- dat1[,!grepl("*bindSpike",names(dat1))]
dat1_pseudoneutid50 <- dat1_pseudoneutid50[,!grepl("*bindRBD",names(dat1_pseudoneutid50))]
dat1_pseudoneutid50 <- dat1_pseudoneutid50[,!grepl("*pseudoneutid80",names(dat1_pseudoneutid50))]



dat1_pseudoneutid50 <- dat1_pseudoneutid50 %>%
  mutate(assay = "pseudoneutid50") %>%
  rename(B = Bpseudoneutid50,
         Day29 = Day29pseudoneutid50,
         Day57 = Day57pseudoneutid50,
         Delta29overB = Delta29overBpseudoneutid50,
         Delta57overB = Delta57overBpseudoneutid50,
         Delta57over29 = Delta57over29pseudoneutid50
  )

dat1_pseudoneutid80 <- dat1[,!grepl("*bindSpike",names(dat1))]
dat1_pseudoneutid80 <- dat1_pseudoneutid80[,!grepl("*bindRBD",names(dat1_pseudoneutid80))]
dat1_pseudoneutid80 <- dat1_pseudoneutid80[,!grepl("*pseudoneutid50",names(dat1_pseudoneutid80))]



dat1_pseudoneutid80 <- dat1_pseudoneutid80 %>%
  mutate(assay = "pseudoneutid80") %>%
  rename(B = Bpseudoneutid80,
         Day29 = Day29pseudoneutid80,
         Day57 = Day57pseudoneutid80,
         Delta29overB = Delta29overBpseudoneutid80,
         Delta57overB = Delta57overBpseudoneutid80,
         Delta57over29 = Delta57over29pseudoneutid80
  )

dat.long <- rbind(dat1_bindSpike,dat1_bindRBD,dat1_pseudoneutid50,dat1_pseudoneutid80) %>% 
  select(Ptid, 
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
         TwophasesampIndD57, 
         ph1.D29,
         ph2.D29,
         Bstratum, 
         wt.D57,  
         wt.D29, 
         wt.intercurrent.cases,
         race, 
         WhiteNonHispanic,
         cohort_event,
         ph1.D57, 
         ph2.D57,
         assay,
         B,
         Day29,
         Day57,
         Delta29overB,
         Delta57overB)


#dat.cor.subset <- dat %>% 
#  filter(TwophasesampIndD57 == 1,
#         ph1.D29 == 1)

#dat.long.cor.subset <- dat.long %>% 
#  filter(TwophasesampIndD57 == 1,
#         ph1.D29 == 1)

dat.long.cor.subset <- dat.long

dat.long.cor.subset <- dat.long.cor.subset %>% 
  mutate(Dich_RaceEthnic = ifelse(EthnicityHispanic == 1, "Hispanic or Latino",
                                  ifelse(EthnicityHispanic == 0 & EthnicityNotreported == 0 & EthnicityUnknown == 0, "Not Hispanic or Latino", NA)
  ))


dat.long.cor.subset <- dat.long.cor.subset %>% mutate(LLoD = ifelse(assay == "pseudoneutid50", log10(2.42),
                                                                    ifelse(assay == "pseudoneutid80", log10(15.02),
                                                                           ifelse(assay == "bindSpike", log10(0.3076),
                                                                                  ifelse(assay == "bindRBD", log10(1.593648), NA)
                                                                           ))))%>% 
                                                mutate(LLoQ = ifelse(assay == "pseudoneutid50", log10(4.477),
                                                                     ifelse(assay == "pseudoneutid80", log10(21.4786),
                                                                            ifelse(assay == "bindSpike", log10(1.7968),
                                                                                   ifelse(assay == "bindRBD", log10(3.4263), NA)
                                                                            ))))%>%
                                                mutate(pos.cutoffs = ifelse(assay == "pseudoneutid50", NA,
                                                                     ifelse(assay == "pseudoneutid80", NA,
                                                                            ifelse(assay == "bindSpike", log10(10.8424),
                                                                                   ifelse(assay == "bindRBD", log10(14.0858), NA)
                                                                            ))))%>%
                                                mutate(ULoQ = ifelse(assay == "pseudoneutid50", log10(10919),
                                                                     ifelse(assay == "pseudoneutid80", log10(15368),
                                                                            ifelse(assay == "bindSpike", log10(10155.95),
                                                                                   ifelse(assay == "bindRBD", log10(16269.23), NA)
                                                                            ))))

dat.long.cor.subset$B[dat.long.cor.subset$B>dat.long.cor.subset$ULoQ] <- dat.long.cor.subset$ULoQ[dat.long.cor.subset$B>dat.long.cor.subset$ULoQ]
dat.long.cor.subset$Day29[dat.long.cor.subset$Day29>dat.long.cor.subset$ULoQ] <- dat.long.cor.subset$ULoQ[dat.long.cor.subset$Day29>dat.long.cor.subset$ULoQ]
dat.long.cor.subset$Day57[dat.long.cor.subset$Day57>dat.long.cor.subset$ULoQ] <- dat.long.cor.subset$ULoQ[dat.long.cor.subset$Day57>dat.long.cor.subset$ULoQ]

dat.long.cor.subset <- dat.long.cor.subset %>% mutate(Delta29overB = Day29 - B, Delta57overB = Day57 - B)




dat.long.cor.subset <- dat.long.cor.subset %>%
  mutate(    
    age_geq_65_label = ifelse(age.geq.65 == 1, "Age >= 65", "Age < 65"),
    highrisk_label = ifelse(HighRiskInd == 1, "At risk", "Not at risk"),
    sex_label = ifelse(Sex == 1, "Female", "Male"),
    ethnicity_label = ifelse(EthnicityHispanic == 1, "Hispanic or Latino", ifelse(EthnicityHispanic == 0 & EthnicityNotreported == 0 & EthnicityUnknown == 0, "Not Hispanic or Latino", "Not reported and unknown")),
    minority_label = ifelse(WhiteNonHispanic == 0, "Comm. of Color", ifelse(WhiteNonHispanic == 1, "White Non-Hispanic", "")),
    Trt = ifelse(Trt == 1,"Vaccine",ifelse(Trt == 0,"Placebo",NA)),
    Bserostatus = ifelse(Bserostatus == 1, "Baseline Pos",ifelse(Bserostatus == 0,"Baseline Neg",NA))
  )

dat.long.cor.subset <- dat.long.cor.subset %>%
  mutate(demo_lab = paste(dat.long.cor.subset$age_geq_65_label, tolower(dat.long.cor.subset$highrisk_label), sep = " ")) %>%
  mutate(demo_lab = as.factor(demo_lab))

dat.long.cor.subset <- dat.long.cor.subset %>%
  mutate(trt_bstatus_label = paste(dat.long.cor.subset$Trt, dat.long.cor.subset$Bserostatus, sep = " ")) %>%
  mutate(trt_bstatus_label = as.factor(trt_bstatus_label))


dat.long.cor.subset <- dat.long.cor.subset %>%
  mutate(age_risk_label = paste(dat.long.cor.subset$age_geq_65_label, tolower(dat.long.cor.subset$highrisk_label), sep = " ")) %>%
  mutate(age_sex_label = paste(dat.long.cor.subset$age_geq_65_label, dat.long.cor.subset$sex_label, sep = " ")) %>%
  mutate(age_minority_label = paste(dat.long.cor.subset$age_geq_65_label, dat.long.cor.subset$minority_label, sep = " ")
  )

dat.long.cor.subset.twophase.intercurrent <- dat.long.cor.subset %>% 
  filter(ph2.D29 == 1) %>%
  filter(!(cohort_event %in% c("Post-Peak Cases","Non-Cases") & ph2.D57==0))

dat.long.cor.subset <- dat.long.cor.subset %>% 
  filter(ph2.D57 == 1)

dat.cor.subset <- dat %>% 
  filter(ph2.D57 == 1)


dat.longer.cor.subset <- dat.long.cor.subset.twophase.intercurrent %>% select(
  Ptid, Trt, Bserostatus, EventIndPrimaryD29, EventIndPrimaryD57, Perprotocol,
  cohort_event, Age, age_geq_65_label, highrisk_label, age_risk_label, sex_label,
  minority_label, Dich_RaceEthnic, assay, LLoD, LLoQ, ULoQ, wt.D57, wt.D29, wt.intercurrent.cases, B, Day29, Day57, Delta29overB, Delta57overB, pos.cutoffs,ph2.D57
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
  mutate(baseline_lt_thres = ifelse(time == "Day 1" & value >= LLoQ, 1, 0)) %>%
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
    response_nab = ifelse((baseline_lt_thres_ptid == 0 & value >= LLoQ) |
                        (baseline_lt_thres_ptid == 1 & time == "Day 1") |
                        (baseline_lt_thres_ptid == 1 & time == "Day 29" & increase_4F_D29_ptid == 1) |
                        (baseline_lt_thres_ptid == 1 & time == "Day 57" & increase_4F_D57_ptid == 1), 1, 0)) %>%
  mutate(
    response_bind = ifelse(value >= pos.cutoffs,1,0)) %>%
  mutate(
    response = case_when(
      assay == "bindSpike" ~ response_bind,
      assay == "bindRBD" ~ response_bind,
      assay == "bindN" ~ response_bind,
      assay == "pseudoneutid50" ~ response_nab,
      assay == "pseudoneutid80" ~ response_nab
    )) %>%
  mutate(
    wt.D29.1 = ifelse(cohort_event == "Intercurrent Cases",wt.intercurrent.cases,wt.D57)
  )




dat.longer.cor.subset.plot1_verification <- dat.longer.cor.subset %>% 
  
  group_by(Trt, Bserostatus, cohort_event, time, assay) %>%
  mutate(num = round(sum(response * ifelse(cohort_event == "Intercurrent Cases", 1, wt.D57)),1),
         denom = round(sum(ifelse(cohort_event=="Intercurrent Cases", 1, wt.D57)),1),
         RespRate = paste(num,"/",denom,"=",round(num/denom*100,1),"%",sep = ""),
         N_RespRate = paste(n(),"\n",round(num/denom*100,1),"%",sep = ""),
         min = min(value),
         q1 = quantile(value, 0.25),
         median = median(value),
         q3 = quantile(value, 0.75),
         max = max(value)) %>%
  ungroup()%>% 
  arrange(Ptid,assay)%>%
  select(-wt.D29.1)



write.csv(dat.longer.cor.subset.plot1_verification,
          here("cor_graphical/verification/output/dat.longer.cor.subset.plot1_verification.csv"),
          row.names = FALSE
)


plot.25sample1<-read.csv(
  here("cor_graphical/verification/input/plot.25sample1.csv")
)

#the subjects plotted in the line plot were randomly selected by original programmer. To verify the input data for line plots
#generated by the original programmer, get the list of combination of subjects and assay that were selected by original programmer.
id_assay <- paste(plot.25sample1$Ptid,plot.25sample1$assay)

plot.25sample1_verification <- dat.longer.cor.subset.plot1_verification[paste(dat.longer.cor.subset.plot1_verification$Ptid,dat.longer.cor.subset.plot1_verification$assay) %in% c(unlist(id_assay, use.names = FALSE)),] %>%
  filter(time %in% c("Day 1","Day 29","Day 57"))

#add "/n" to the response rate to match the print out format shown in the plot.
plot.25sample1_verification$RespRate <- gsub("=", "\n", plot.25sample1_verification$RespRate )


write.csv(plot.25sample1_verification,
          here("cor_graphical/verification/output/plot.25sample1_verification.csv"),
          row.names = FALSE
)


#lineplots of Binding Antibody to Spike: baseline negative vaccine arm (2 timepoints)
lineplots_neg_vaccine_bindSpike_2 <- dat.longer.cor.subset.plot1_verification %>% 
  filter(Trt == "Vaccine", 
         Bserostatus == "Baseline Neg", 
         assay == "bindSpike",
         time %in% c("Day 29","Day 57"))

#the subjects plotted in the line plot were randomly selected by original programmer. To verify the line plots
#generated by the original programmer, get the list of combination of subjects and assay that were selected by original programmer.
id <-plot.25sample1 %>%
  filter(Trt == "Vaccine", 
         Bserostatus == "Baseline Neg", 
         assay == "bindSpike",
         time %in% c("Day 29 ","Day 57")) %>%
  select(Ptid)

lineplots_neg_vaccine_bindSpike_plot.25sample1_2 <- lineplots_neg_vaccine_bindSpike_2[lineplots_neg_vaccine_bindSpike_2$Ptid %in% c(unlist(id, use.names = FALSE)),]

p <- ggplot(lineplots_neg_vaccine_bindSpike_2, aes(x = time, y = value)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1,coef = 6) +
  scale_y_continuous(breaks = c(-1,0, 1, 2, 3, 4, 5, 6, 7, 8, 9), limits = c(-1, 10)) +
  facet_wrap(~cohort_event) +
  geom_line(data=lineplots_neg_vaccine_bindSpike_plot.25sample1_2,aes(x = time, y = value, group=Ptid)) +
  geom_hline(yintercept=lineplots_neg_vaccine_bindSpike_2$pos.cutoffs, linetype="dashed", color = "black") +
  geom_hline(yintercept=lineplots_neg_vaccine_bindSpike_2$ULoQ, linetype="dashed", color = "black")

p

lineplots_neg_vaccine_bindSpike_2 %>% 
  select(cohort_event, time, RespRate, N_RespRate, min, q1, median, q3, max) %>% 
  unique()


#lineplots of Pseudovirus Neutralization ID50: baseline negative vaccine arm (2 timepoints)
lineplots_neg_vaccine_pseudoneutid50_2 <- dat.longer.cor.subset.plot1_verification %>% 
  filter(Trt == "Vaccine", 
         Bserostatus == "Baseline Neg", 
         assay == "pseudoneutid50",
         time %in% c("Day 29","Day 57"))

#the subjects plotted in the line plot were randomly selected by original programmer. To verify the line plots
#generated by the original programmer, get the list of combination of subjects and assay that were selected by original programmer.
id <-plot.25sample1 %>%
  filter(Trt == "Vaccine", 
         Bserostatus == "Baseline Neg", 
         assay == "pseudoneutid50",
         time %in% c("Day 29","Day 57")) %>%
  select(Ptid)

lineplots_neg_vaccine_pseudoneutid50_plot.25sample1_2 <- lineplots_neg_vaccine_pseudoneutid50_2[lineplots_neg_vaccine_pseudoneutid50_2$Ptid %in% c(unlist(id, use.names = FALSE)),]

p <- ggplot(lineplots_neg_vaccine_pseudoneutid50_2, aes(x = time, y = value)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1,coef = 6) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7), limits = c(0, 10)) +
  facet_wrap(~cohort_event) +
  geom_line(data=lineplots_neg_vaccine_pseudoneutid50_plot.25sample1_2,aes(x = time, y = value, group=Ptid)) +
  geom_hline(yintercept=lineplots_neg_vaccine_pseudoneutid50_2$LLoD, linetype="dashed", color = "black") +
  geom_hline(yintercept=lineplots_neg_vaccine_pseudoneutid50_2$ULoQ, linetype="dashed", color = "black")

p

lineplots_neg_vaccine_pseudoneutid50_2 %>% 
  select(cohort_event, time, RespRate, N_RespRate, min, q1, median, q3, max) %>% 
  unique()


#lineplots of Binding Antibody to Spike: baseline negative vaccine arm (3 timepoints)
lineplots_neg_vaccine_bindSpike_3 <- dat.longer.cor.subset.plot1_verification %>% 
  filter(Trt == "Vaccine", 
         Bserostatus == "Baseline Neg", 
         assay == "bindSpike",
         time %in% c("Day 1","Day 29","Day 57"))

#the subjects plotted in the line plot were randomly selected by original programmer. To verify the line plots
#generated by the original programmer, get the list of combination of subjects and assay that were selected by original programmer.
id <-plot.25sample1 %>%
  filter(Trt == "Vaccine", 
         Bserostatus == "Baseline Neg", 
         assay == "bindSpike",
         time %in% c("Day 1","Day 29","Day 57")) %>%
  select(Ptid)


lineplots_neg_vaccine_bindSpike_plot.25sample1_3 <- lineplots_neg_vaccine_bindSpike_3[lineplots_neg_vaccine_bindSpike_3$Ptid %in% c(unlist(id, use.names = FALSE)),]

p <- ggplot(lineplots_neg_vaccine_bindSpike_3, aes(x = time, y = value)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1,coef = 6) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7), limits = c(-1, 10)) +
  facet_wrap(~cohort_event) +
  geom_line(data=lineplots_neg_vaccine_bindSpike_plot.25sample1_3,aes(x = time, y = value, group=Ptid)) +
  geom_hline(yintercept=lineplots_neg_vaccine_bindSpike_3$pos.cutoffs, linetype="dashed", color = "black") +
  geom_hline(yintercept=lineplots_neg_vaccine_bindSpike_3$ULoQ, linetype="dashed", color = "black")

p

lineplots_neg_vaccine_bindSpike_3 %>% 
  select(cohort_event, time, RespRate, N_RespRate, min, q1, median, q3, max) %>% 
  unique()




