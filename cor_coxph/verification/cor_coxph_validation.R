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
inputFile <- read.csv(here::here("..", "data_clean", data_name)) 

# Define marker data for analysis, Day 57 or Day 29 as this will determine the data subset for analysis (Add the intercurrent cases in Day 29 analyses!)
markerDay = 57

# create data subsets
phase1.D57 <- inputFile %>% filter(is.na(wt))
phase2.D57 <- inputFile %>% filter(TwophasesampInd == 1)
dat.D57 <- bind_rows(phase1.D57, phase2.D57)

phase1.D29 <- inputFile %>% filter(is.na(wt.2))
phase2.D29 <- inputFile %>% filter(TwophasesampInd.2 == 1)
dat.D29 <- bind_rows(phase1.D29, phase2.D29)

# Define trichotomized version of the markers: Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80
vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D57$Day57bindSpike, probs = c(1/3, 2/3)) + c(1e-6, 0), Inf)
dat.D57$Day57bindSpikecat <- Hmisc::cut2(dat.D57$Day57bindSpike, vec_cutpoints)
rm(vec_cutpoints)

vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D57$Day57bindRBD, probs = c(1/3, 2/3)) + c(1e-6, 0), Inf)
dat.D57$Day57bindRBDcat <- Hmisc::cut2(dat.D57$Day57bindRBD, vec_cutpoints)
rm(vec_cutpoints)

vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D57$Day57pseudoneutid50, probs = c(1/3, 2/3)) + c(1e-6, 0), Inf)
dat.D57$Day57pseudoneutid50cat <- Hmisc::cut2(dat.D57$Day57pseudoneutid50, vec_cutpoints)
rm(vec_cutpoints)

vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D57$Day57pseudoneutid80, probs = c(1/3, 2/3)) + c(1e-6, 0), Inf)
dat.D57$Day57pseudoneutid80cat <- Hmisc::cut2(dat.D57$Day57pseudoneutid80, vec_cutpoints)
rm(vec_cutpoints)

dat.D57.design <- twophase(id=list(~Ptid,~Ptid),strata=list(NULL,~Wstratum),
                           data=dat.D57, subset=~TwophasesampInd)

fit <- survey::svycoxph(Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ Day57bindSpike + Day57bindRBD + Day57pseudoneutid50 + Day57pseudoneutid80 + MinorityInd + HighRiskInd + Age, design=dat.D57.design)

survey::svycoxph(Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ Day57bindSpikecat + Day57bindRBDcat + Day57pseudoneutid50cat + Day57pseudoneutid80cat + MinorityInd + HighRiskInd + Age, design=dat.D57.design)


