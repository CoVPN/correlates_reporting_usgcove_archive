###########################
# File: immuno_graphical_verification_code.R
# Author: Di Lu
# Creation Date: 2021/03/19
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
library(here)


# initial data processing
dat <- read.csv(
  here("immuno_graphical/verification/verification_input/practice_data.csv")
)


data.variables_needed <- dat %>% 
  select(
  Ptid, 
  Trt, 
  MinorityInd, 
  HighRiskInd, 
  Age, 
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
  wt.subcohort, 
  race,
  ethnicity, 
  EthnicityHispanic, 
  EthnicityNotreported, 
  EthnicityUnknown, 
  WhiteNonHispanic,
  BbindSpike, 
  BbindRBD, 
  Bpseudoneutid50, 
  Bpseudoneutid80,
  Day29bindSpike, 
  Day29bindRBD, 
  Day29pseudoneutid50, 
  Day29pseudoneutid80,
  Day57bindSpike, 
  Day57bindRBD, 
  Day57pseudoneutid50,
  Day57pseudoneutid80,
  Delta29overBbindSpike, 
  Delta29overBbindRBD, 
  Delta29overBpseudoneutid50, 
  Delta29overBpseudoneutid80,
  Delta57overBbindSpike, 
  Delta57overBbindRBD, 
  Delta57overBpseudoneutid50, 
  Delta57overBpseudoneutid80,
  Delta57over29bindSpike, 
  Delta57over29bindRBD, 
  Delta57over29pseudoneutid50, 
  Delta57over29pseudoneutid80
)

# wide to long
data.variables_needed_bindSpike <- dat %>%
  select(
    Ptid, 
    Trt, 
    MinorityInd, 
    HighRiskInd, 
    Age, 
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
    wt.subcohort, 
    race, 
    ethnicity, 
    EthnicityHispanic, 
    EthnicityNotreported,
    EthnicityUnknown, 
    WhiteNonHispanic,
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
    HighRiskInd, 
    Age, 
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
    wt.subcohort, 
    race, 
    ethnicity, 
    EthnicityHispanic, 
    EthnicityNotreported,
    EthnicityUnknown, 
    WhiteNonHispanic,
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
    HighRiskInd, 
    Age, 
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
    wt.subcohort, 
    race, 
    ethnicity, 
    EthnicityHispanic, 
    EthnicityNotreported,
    EthnicityUnknown, 
    WhiteNonHispanic,
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
    HighRiskInd, 
    Age, 
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
    wt.subcohort, 
    race, 
    ethnicity, 
    EthnicityHispanic, 
    EthnicityNotreported,
    EthnicityUnknown, 
    WhiteNonHispanic,
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



dat1 <- rbind(
  data.variables_needed_bindSpike, 
  data.variables_needed_bindRBD, 
  data.variables_needed_pseudoneutid50, 
  data.variables_needed_pseudoneutid80)

# keep the cohort identified for immunogenicity
dat.long.twophase.sample <- dat1 %>% 
  filter(SubcohortInd == 1, 
         TwophasesampInd == 1, 
         Perprotocol == 1)

dat.long.twophase.sample <- dat.long.twophase.sample %>%
  mutate(
    trt_bstatus_label = paste(dat.long.twophase.sample$Trt, dat.long.twophase.sample$Bserostatus, sep = "*"),
    age_geq_65_label = ifelse(age.geq.65 == 1, ">=65", "<65"),
    highrisk_label = ifelse(HighRiskInd == 1, "at risk", "not at risk"),
    sex_label = ifelse(Sex == 1, "Male", "Female"),
    ethnicity_label = ifelse(EthnicityHispanic == 1, "Hispanic or Latino", ifelse(EthnicityHispanic == 0 & EthnicityNotreported == 0 & EthnicityUnknown == 0, "Not Hispanic or Latino", "Not reported and unknown")),
    minority_label = ifelse(WhiteNonHispanic == 0, "Comm. of Color", ifelse(WhiteNonHispanic == 1, "White Non-Hispanic", ""))
    )
 

dat.long.twophase.sample <- dat.long.twophase.sample %>%
  mutate(age_risk_label = paste(dat.long.twophase.sample$age_geq_65_label, dat.long.twophase.sample$highrisk_label, sep = "*"),
         age_sex_label = paste(dat.long.twophase.sample$age_geq_65_label, dat.long.twophase.sample$sex_label, sep = "*"),
         age_minority_label = paste(dat.long.twophase.sample$age_geq_65_label, dat.long.twophase.sample$minority_label, sep = "*")
         ) 

#Pair plots of D57 fold-rise over D1 Ab markers: baseline negative vaccine arm
bnegative_vaccine_long <- dat.long.twophase.sample %>%
  filter(Trt == 1, Bserostatus == 0) %>%
  select(Ptid, assay, Day57)

bnegative_vaccine_wide <- reshape(bnegative_vaccine_long,
  timevar = "assay",
  idvar = "Ptid",
  direction = "wide"
)


pairs.panels(bnegative_vaccine_wide[, c(2:5)],
  method = "spearman",
  hist.col = "#00AFBB",
  density = TRUE,
  ellipses = TRUE,
  xlim = c(0, 10),
  ylim = c(0, 10)
)

# figure 1.3 Verify the Partial Spearman's Rank Correlation
spearman_corr <- dat.long.twophase.sample %>%
  filter(
    Trt == 1, 
    Bserostatus == 0) %>%
  select(Ptid, 
         assay, 
         Day57, 
         Bstratum, 
         wt.subcohort)

spearman_corr_wide <- reshape(spearman_corr,
  timevar = "assay",
  idvar = c("Ptid", "Bstratum", "wt.subcohort"),
  direction = "wide") %>%
  mutate(x1 = ifelse(Bstratum == 2, 1, 0)) %>%
  mutate(x2 = ifelse(Bstratum == 3, 1, 0))

spearman_corr_wide <- spearman_corr_wide %>% 
  drop_na(wt.subcohort)

spearman_corr_wide$wt.subcohort <- spearman_corr_wide$wt.subcohort / sum(spearman_corr_wide$wt.subcohort)

spearman_corr_calculation <- function(var1_var2) {
  c <- 0
  for (i in 1:500) {
    resample_id <- sample(x = c(1:nrow(spearman_corr_wide)), size = nrow(spearman_corr_wide), replace = TRUE, prob = spearman_corr_wide$wt.subcohort)
    resample <- spearman_corr_wide[resample_id, ]
    formula_x <- formula(paste0(var1_var2, "~", paste0("x", 1:2, collapse = "+")))
    a <- partial_Spearman(formula_x, fit.x = "lm", fit.y = "lm", data = resample)
    c[i] <- a$TS$TB$ts
  }
  return(list(var1_var2, mean(c)))
}

set.seed(12345)
spearman_corr_calculation("Day57.bindSpike|Day57.bindRBD")
set.seed(12345)
spearman_corr_calculation("Day57.bindSpike|Day57.pseudoneutid50")
set.seed(12345)
spearman_corr_calculation("Day57.bindSpike|Day57.pseudoneutid80")
set.seed(12345)
spearman_corr_calculation("Day57.bindRBD|Day57.pseudoneutid50")
set.seed(12345)
spearman_corr_calculation("Day57.bindRBD|Day57.pseudoneutid80")
set.seed(12345)
spearman_corr_calculation("Day57.pseudoneutid50|Day57.pseudoneutid80")


#RCDF plots for D57 Ab markers: baseline negative vaccine arm
RCDF_plots_data_neg_long <- dat.long.twophase.sample %>%
  filter(Trt == 1, 
         Bserostatus == 0) %>%
  select(Ptid, 
         assay, 
         Day57, 
         wt.2) %>%
  drop_na(wt.2)

RCDF_plots_data_pos_long <- dat.long.twophase.sample %>%
  filter(Trt == 1,
         Bserostatus == 1) %>%
  select(Ptid, 
         assay, 
         Day57, 
         wt.2) %>%
  drop_na(wt.2)

RCDF <- function(assay_type) {
  df <- RCDF_plots_data %>% 
    filter(assay == assay_type)
  df <- df[order(df$Day57), ]
  df$cum.pct <- with(df, cumsum(wt.2) / sum(wt.2))
  df$r_cum.pct <- 1 - df$cum.pct

  start <- floor(df[1, ]$Day57)
  end <- ceiling(df[nrow(df), ]$Day57)
  if (start > 0) {
    for (i in 1:start) {
      print(1)
    }
  }
  for (j in (start + 1):(end - 1)) {
    print(df[which.min(abs(df$Day57 - j)), ]$r_cum.pct)
  }
  for (k in end:10) {
    print(0)
  }
}

RCDF_plots_data <- RCDF_plots_data_neg_long
RCDF("bindSpike")
RCDF("bindRBD")
RCDF("pseudoneutid50")
RCDF("pseudoneutid80")

RCDF_plots_data <- RCDF_plots_data_pos_long
RCDF("bindSpike")
RCDF("bindRBD")
RCDF("pseudoneutid50")
RCDF("pseudoneutid80")

RCDF_step <- dat.long.twophase.sample %>% 
  filter(Trt == 1) %>%
  select(Ptid, 
         assay, 
         Bserostatus,
         Day57, 
         wt.2) %>%
  mutate(assay_Bserostatus = paste(assay,Bserostatus, sep = "*"))


RCDF_step_bindSpike_1 <- RCDF_step %>% filter(assay == "bindSpike",Bserostatus == 1)
RCDF_step_bindSpike_1 <- RCDF_step_bindSpike_1[order(RCDF_step_bindSpike_1$Day57), ]
RCDF_step_bindSpike_1$cum.pct <- with(RCDF_step_bindSpike_1, 1-cumsum(wt.2) / sum(wt.2))

RCDF_step_bindSpike_0 <- RCDF_step %>% filter(assay == "bindSpike",Bserostatus == 0)
RCDF_step_bindSpike_0 <- RCDF_step_bindSpike_0[order(RCDF_step_bindSpike_0$Day57), ]
RCDF_step_bindSpike_0$cum.pct <- with(RCDF_step_bindSpike_0, 1-cumsum(wt.2) / sum(wt.2))

RCDF_step_bindRBD_1<-RCDF_step%>%filter(assay == "bindRBD",Bserostatus == 1)
RCDF_step_bindRBD_1 <- RCDF_step_bindRBD_1[order(RCDF_step_bindRBD_1$Day57), ]
RCDF_step_bindRBD_1$cum.pct <- with(RCDF_step_bindRBD_1, 1-cumsum(wt.2) / sum(wt.2))

RCDF_step_bindRBD_0<-RCDF_step%>%filter(assay == "bindRBD",Bserostatus == 0)
RCDF_step_bindRBD_0 <- RCDF_step_bindRBD_0[order(RCDF_step_bindRBD_0$Day57), ]
RCDF_step_bindRBD_0$cum.pct <- with(RCDF_step_bindRBD_0, 1-cumsum(wt.2) / sum(wt.2))

RCDF_step_pseudoneutid50_1<-RCDF_step%>%filter(assay == "pseudoneutid50",Bserostatus == 1)
RCDF_step_pseudoneutid50_1 <- RCDF_step_pseudoneutid50_1[order(RCDF_step_pseudoneutid50_1$Day57), ]
RCDF_step_pseudoneutid50_1$cum.pct <- with(RCDF_step_pseudoneutid50_1, 1-cumsum(wt.2) / sum(wt.2))

RCDF_step_pseudoneutid50_0<-RCDF_step%>%filter(assay == "pseudoneutid50",Bserostatus == 0)
RCDF_step_pseudoneutid50_0 <- RCDF_step_pseudoneutid50_0[order(RCDF_step_pseudoneutid50_0$Day57), ]
RCDF_step_pseudoneutid50_0$cum.pct <- with(RCDF_step_pseudoneutid50_0, 1-cumsum(wt.2) / sum(wt.2))

RCDF_step_pseudoneutid80_1<-RCDF_step%>%filter(assay == "pseudoneutid80",Bserostatus == 1)
RCDF_step_pseudoneutid80_1 <- RCDF_step_pseudoneutid80_1[order(RCDF_step_pseudoneutid80_1$Day57), ]
RCDF_step_pseudoneutid80_1$cum.pct <- with(RCDF_step_pseudoneutid80_1, 1-cumsum(wt.2) / sum(wt.2))

RCDF_step_pseudoneutid80_0<-RCDF_step%>%filter(assay == "pseudoneutid80",Bserostatus == 0)
RCDF_step_pseudoneutid80_0 <- RCDF_step_pseudoneutid80_0[order(RCDF_step_pseudoneutid80_0$Day57), ]
RCDF_step_pseudoneutid80_0$cum.pct <- with(RCDF_step_pseudoneutid80_0, 1-cumsum(wt.2) / sum(wt.2))

RCDF_step<-rbind(RCDF_step_bindSpike_1,
                 RCDF_step_bindSpike_0,
                 RCDF_step_bindRBD_1,
                 RCDF_step_bindRBD_0,
                 RCDF_step_pseudoneutid50_1,
                 RCDF_step_pseudoneutid50_0,
                 RCDF_step_pseudoneutid80_1,
                 RCDF_step_pseudoneutid80_0)
        
ggplot(RCDF_step, aes(Day57, cum.pct,group=assay_Bserostatus, color = assay_Bserostatus),xlim=c(-2,10)) + 
  geom_step() +
  scale_x_continuous(breaks = c(0,2,4,6,8,10))+
  ylab("Reverse ECDF")

#Boxplots of D29 Ab markers: baseline negative vaccine + placebo arms

boxplots_data_long <- dat.long.twophase.sample %>% select(Ptid, Trt, Bserostatus, assay, Day57)
# baseline negtive placebo arm
boxplots_data_long %>%
  filter(
    Trt == 0, 
    Bserostatus == 0) %>%
  group_by(assay) %>%
  summarize(
    min = min(Day57),
    q1 = quantile(Day57, 0.25),
    median = quantile(Day57, 0.50),
    mean = mean(Day57),
    q3 = quantile(Day57, 0.75),
    max = max(Day57),
    miniqr = quantile(Day57, 0.25) - 1.5 * (quantile(Day57, 0.75) - quantile(Day57, 0.25)),
    maxiqr = quantile(Day57, 0.75) + 1.5 * (quantile(Day57, 0.75) - quantile(Day57, 0.25))
  )
# baseline Negative vaccine arm
boxplots_data_long %>%
  filter(
    Trt == 1, 
    Bserostatus == 0) %>%
  group_by(assay) %>%
  summarize(
    min = min(Day57),
    q1 = quantile(Day57, 0.25),
    median = quantile(Day57, 0.50),
    mean = mean(Day57),
    q3 = quantile(Day57, 0.75),
    max = max(Day57),
    miniqr = quantile(Day57, 0.25) - 1.5 * (quantile(Day57, 0.75) - quantile(Day57, 0.25)),
    maxiqr = quantile(Day57, 0.75) + 1.5 * (quantile(Day57, 0.75) - quantile(Day57, 0.25))
  ) 

# baseline Negative placebo arm
boxplots_data_long1<-boxplots_data_long
boxplots_data_long1$Bserostatus<-as.factor(boxplots_data_long1$Bserostatus)
boxplots_data_long1<-boxplots_data_long1 %>% 
  filter(Bserostatus==0,
         Trt == 0)

ggplot(boxplots_data_long1, aes(y = Day57)) +
  stat_boxplot(aes(x=assay),
               geom = "errorbar", 
               width = 0.5) +  
  geom_boxplot(aes(x=assay))  +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10))

# baseline Negative vaccine arm
boxplots_data_long1<-boxplots_data_long
boxplots_data_long1$Bserostatus<-as.factor(boxplots_data_long1$Bserostatus)
boxplots_data_long1<-boxplots_data_long1 %>% 
  filter(Bserostatus==0,
         Trt == 1)

ggplot(boxplots_data_long1, aes(y = Day57)) +
  stat_boxplot(aes(x=assay),
               geom = "errorbar", 
               width = 0.5) +  
  geom_boxplot(aes(x=assay))  +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10))



#Spaghetti plots of Ab markers over time: baseline negative vaccine + placebo arm

Spaghetti_data_long <- dat.long.twophase.sample %>% 
  select(Ptid, 
         Trt, 
         Bserostatus, 
         assay, 
         B, 
         Day29, 
         Day57)

Spaghetti_data_long %>%
  filter(
    Trt == 0, 
    Bserostatus == 0) %>%
  group_by(assay) %>%
  summarize(
    min_B = min(B),
    q1_B = quantile(B, 0.25),
    median_B = quantile(B, 0.50),
    mean_B = mean(B),
    q3_B = quantile(B, 0.75),
    max_B = max(B),
    
    min_D29 = min(Day29),
    q1_D29 = quantile(Day29, 0.25),
    median_D29 = quantile(Day29, 0.50),
    mean_D29 = mean(Day29),
    q3_D29 = quantile(Day29, 0.75),
    max_D29 = max(Day29),
    
    min_D57 = min(Day57),
    q1_D57 = quantile(Day57, 0.25),
    median_D57 = quantile(Day57, 0.50),
    mean_D57 = mean(Day57),
    q3_D57 = quantile(Day57, 0.75),
    max_D57 = max(Day57)
  )





Spaghetti_data_long <- dat.long.twophase.sample %>% 
  select(
    Ptid, 
    Trt, 
    Bserostatus, 
    assay, 
    B, 
    Day29, 
    Day57)

Spaghetti_data_long %>%
  filter(
    Trt == 1, 
    Bserostatus == 0) %>%
  group_by(assay) %>%
  summarize(
    min_B = min(B),
    q1_B = quantile(B, 0.25),
    median_B = quantile(B, 0.50),
    mean_B = mean(B),
    q3_B = quantile(B, 0.75),
    max_B = max(B),
    
    min_D29 = min(Day29),
    q1_D29 = quantile(Day29, 0.25),
    median_D29 = quantile(Day29, 0.50),
    mean_D29 = mean(Day29),
    q3_D29 = quantile(Day29, 0.75),
    max_D29 = max(Day29),
    
    min_D57 = min(Day57),
    q1_D57 = quantile(Day57, 0.25),
    median_D57 = quantile(Day57, 0.50),
    mean_D57 = mean(Day57),
    q3_D57 = quantile(Day57, 0.75),
    max_D57 = max(Day57)
  )



