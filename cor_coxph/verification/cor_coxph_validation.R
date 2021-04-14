#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

# load required libraries and functions
library(tidyverse)
library(here)
library(methods)
library(SuperLearner)
library(glmnet)
library(kyotil)
library(conflicted)
conflicted::conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
library(mice)
library(marginalizedRisk)
library(tools) 
library(survey)
library(parallel)
library(forestplot)
library(Hmisc)
library(xtable) 

# Read in data file
inputFile <- read.csv(here::here("..", "data_clean", "practice_data_with_riskscore.csv")) 

# Define marker data for analysis, Day 57 or Day 29 as this will determine the data subset for analysis (Add the intercurrent cases in Day 29 analyses!)
markerDay = 57

# create data subsets
dat <- inputFile %>% filter(Trt == 1 & Bserostatus == 0)


################################################## 
# Day 57 analysis

dat.D57 <- dat %>% filter(!is.na(wt))

dat.D57.design <- twophase(id=list(~Ptid, ~Ptid),strata=list(NULL, ~Wstratum),
                           data=dat.D57, subset=~TwophasesampInd)

fit <- survey::svycoxph(Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ Day57bindSpike + MinorityInd + HighRiskInd + risk_score, design=dat.D57.design)

# Define trichotomized version of the markers: Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80
vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D57$Day57bindSpike, weights = dat.D57$wt, probs = c(1/3, 2/3)), Inf)
if((vec_cutpoints[[2]] == min(dat.D57$Day57bindSpike, na.rm = T)) | 
   (vec_cutpoints[[3]] == max(dat.D57$Day57bindSpike, na.rm = T)) |
   (vec_cutpoints[[2]] == vec_cutpoints[[3]])){
  dat.D57$Day57bindSpikecat <- Hmisc::cut2(dat.D57$Day57bindSpike, breaks = 3)
} else {
    dat.D57$Day57bindSpikecat <- Hmisc::cut2(dat.D57$Day57bindSpike, vec_cutpoints)}
rm(vec_cutpoints)



vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D57$Day57bindRBD, weights = dat.D57$wt, probs = c(1/3, 2/3)), Inf)
if((vec_cutpoints[[2]] == min(dat.D57$Day57bindRBD, na.rm = T)) | 
   (vec_cutpoints[[3]] == max(dat.D57$Day57bindRBD, na.rm = T)) |
   (vec_cutpoints[[2]] == vec_cutpoints[[3]])){
  dat.D57$Day57bindRBDcat <- Hmisc::cut2(dat.D57$Day57bindRBD, breaks = 3)
} else {
  dat.D57$Day57bindRBDcat <- Hmisc::cut2(dat.D57$Day57bindRBD, vec_cutpoints)}
rm(vec_cutpoints)




vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D57$Day57pseudoneutid50, weights = dat.D57$wt, probs = c(1/3, 2/3)), Inf)
if((vec_cutpoints[[2]] == min(dat.D57$Day57pseudoneutid50, na.rm = T)) | 
   (vec_cutpoints[[3]] == max(dat.D57$Day57pseudoneutid50, na.rm = T)) |
   (vec_cutpoints[[2]] == vec_cutpoints[[3]])){
  dat.D57$Day57pseudoneutid50cat <- Hmisc::cut2(dat.D57$Day57pseudoneutid50, breaks = 3)
} else {
  dat.D57$Day57pseudoneutid50cat <- Hmisc::cut2(dat.D57$Day57pseudoneutid50, vec_cutpoints)}
rm(vec_cutpoints)




vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D57$Day57pseudoneutid80, weights = dat.D57$wt, probs = c(1/3, 2/3)), Inf)
if((vec_cutpoints[[2]] == min(dat.D57$Day57pseudoneutid80, na.rm = T)) | 
   (vec_cutpoints[[3]] == max(dat.D57$Day57pseudoneutid80, na.rm = T)) |
   (vec_cutpoints[[2]] == vec_cutpoints[[3]])){
  dat.D57$Day57pseudoneutid80cat <- cut(dat.D57$Day57pseudoneutid80, breaks = 3)
} else {
  dat.D57$Day57pseudoneutid80cat <- Hmisc::cut2(dat.D57$Day57pseudoneutid80, vec_cutpoints)}
rm(vec_cutpoints)



dat.D57 <- dat.D57 %>%
  mutate(Day57bindSpikecat = factor(Day57bindSpikecat),
         Day57bindRBDcat = factor(Day57bindRBDcat),
         Day57pseudoneutid50cat = factor(Day57pseudoneutid50cat),
         Day57pseudoneutid80cat = factor(Day57pseudoneutid80cat))


dat.D57.design <- twophase(id=list(~Ptid, ~Ptid),strata=list(NULL, ~Wstratum),
                           data=dat.D57, subset=~TwophasesampInd)

fitcat <- survey::svycoxph(Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ Day57bindSpikecat + MinorityInd + HighRiskInd + risk_score, design=dat.D57.design)

################################################## 
# Day 29 analysis
dat.D29 <- dat %>% filter(!is.na(wt.2))

dat.D29.design <- twophase(id=list(~Ptid, ~Ptid),strata=list(NULL, ~Wstratum),
                           data=dat.D29, subset=~TwophasesampInd.2)

fit <- survey::svycoxph(Surv(EventTimePrimaryD29, EventIndPrimaryD29) ~ Day29bindSpike + MinorityInd + HighRiskInd + risk_score, design=dat.D29.design)

# Define trichotomized version of the markers: Day29bindSpike, Day29bindRBD, Day29pseudoneutid50, Day29pseudoneutid80
vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D29$Day29bindSpike, weights = dat.D29$wt.2, probs = c(1/3, 2/3)), Inf)
if((vec_cutpoints[[2]] == min(dat.D29$Day29bindSpike, na.rm = T)) | 
   (vec_cutpoints[[3]] == max(dat.D29$Day29bindSpike, na.rm = T)) |
   (vec_cutpoints[[2]] == vec_cutpoints[[3]])){
  dat.D29$Day29bindSpikecat <- Hmisc::cut2(dat.D29$Day29bindSpike, breaks = 3)
} else {
  dat.D29$Day29bindSpikecat <- Hmisc::cut2(dat.D29$Day29bindSpike, vec_cutpoints)}
rm(vec_cutpoints)



vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D29$Day29bindRBD, weights = dat.D29$wt.2, probs = c(1/3, 2/3)), Inf)
if((vec_cutpoints[[2]] == min(dat.D29$Day29bindRBD, na.rm = T)) | 
   (vec_cutpoints[[3]] == max(dat.D29$Day29bindRBD, na.rm = T)) |
   (vec_cutpoints[[2]] == vec_cutpoints[[3]])){
  dat.D29$Day29bindRBDcat <- Hmisc::cut2(dat.D29$Day29bindRBD, breaks = 3)
} else {
  dat.D29$Day29bindRBDcat <- Hmisc::cut2(dat.D29$Day29bindRBD, vec_cutpoints)}
rm(vec_cutpoints)




vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D29$Day29pseudoneutid50, weights = dat.D29$wt.2, probs = c(1/3, 2/3)), Inf)
if((vec_cutpoints[[2]] == min(dat.D29$Day29pseudoneutid50, na.rm = T)) | 
   (vec_cutpoints[[3]] == max(dat.D29$Day29pseudoneutid50, na.rm = T)) |
   (vec_cutpoints[[2]] == vec_cutpoints[[3]])){
  dat.D29$Day29pseudoneutid50cat <- Hmisc::cut2(dat.D29$Day29pseudoneutid50, breaks = 3)
} else {
  dat.D29$Day29pseudoneutid50cat <- Hmisc::cut2(dat.D29$Day29pseudoneutid50, vec_cutpoints)}
rm(vec_cutpoints)




vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D29$Day29pseudoneutid80, weights = dat.D29$wt.2, probs = c(1/3, 2/3)), Inf)
if((vec_cutpoints[[2]] == min(dat.D29$Day29pseudoneutid80, na.rm = T)) | 
   (vec_cutpoints[[3]] == max(dat.D29$Day29pseudoneutid80, na.rm = T)) |
   (vec_cutpoints[[2]] == vec_cutpoints[[3]])){
  dat.D29$Day29pseudoneutid80cat <- cut(dat.D29$Day29pseudoneutid80, breaks = 3)
} else {
  dat.D29$Day29pseudoneutid80cat <- Hmisc::cut2(dat.D29$Day29pseudoneutid80, vec_cutpoints)}
rm(vec_cutpoints)



dat.D29 <- dat.D29 %>%
  mutate(Day29bindSpikecat = factor(Day29bindSpikecat),
         Day29bindRBDcat = factor(Day29bindRBDcat),
         Day29pseudoneutid50cat = factor(Day29pseudoneutid50cat),
         Day29pseudoneutid80cat = factor(Day29pseudoneutid80cat))


dat.D29.design <- twophase(id=list(~Ptid, ~Ptid),strata=list(NULL, ~Wstratum),
                           data=dat.D29, subset=~TwophasesampInd.2)

fitcat <- survey::svycoxph(Surv(EventTimePrimaryD29, EventIndPrimaryD29) ~ Day29bindSpikecat + MinorityInd + HighRiskInd + risk_score, design=dat.D29.design)


