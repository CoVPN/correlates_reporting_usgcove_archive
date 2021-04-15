#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

# load required libraries and functions
library(tidyverse)
library(here)
# library(methods)
# library(SuperLearner)
library(glmnet)
library(kyotil)
library(conflicted)
conflicted::conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
# library(mice)
library(marginalizedRisk)
library(tools) 
library(survey)
# library(parallel)
# library(forestplot)
library(Hmisc)
library(xtable) 
library(broom)

# Read in data file
inputFile <- read.csv(here::here("..", "data_clean", "practice_data_with_riskscore.csv")) 

# Define marker data for analysis, Day 57 or Day 29 as this will determine the data subset for analysis (Add the intercurrent cases in Day 29 analyses!)
markerDay = 57

# create data subsets
dat <- inputFile %>% filter(Trt == 1 & Bserostatus == 0)

# This function takes marker name as string, and the design.
# It returns the HR for the marker, CI and p.value
getHR_D57_continuous_marker <- function(marker, design, data, group){
  fm = as.formula(paste("Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ ", marker, " + MinorityInd + HighRiskInd + risk_score"))
  
  fit <- survey::svycoxph(fm, design=design)
  
  tidy(fit, conf.int = T, exponentiate = T) %>% data.frame() %>%
    select(term, estimate, conf.low, conf.high, p.value) %>%
    filter(!term %in% c("MinorityInd", "HighRiskInd", "risk_score")) %>%
    rename(marker = term) %>%
    mutate(events = sum(data$EventIndPrimaryD57),
           n = length(data$EventIndPrimaryD57),
           group = group) 
}



getHR_D57_categorical_marker <- function(marker, design, data, group){
  
  # table(data %>% pull(all_of(marker)))
  # Add events and n for lower, middle and upper
  data_lower <- data %>% filter(get(marker) == (data %>% pull(marker) %>% levels())[1])
  data_middle <- data %>% filter(get(marker) == (data %>% pull(marker) %>% levels())[2])
  data_upper <- data %>% filter(get(marker) == (data %>% pull(marker) %>% levels())[3])
  # data_lower <- data %>% filter(str_detect(get(marker), '-Inf'))
  # data_middle <- data %>% filter(!str_detect(get(marker), 'Inf'))
  # data_upper <- data %>% filter(str_detect(get(marker), ', Inf]') | str_detect(marker, '3.11'))

  fm = as.formula(paste("Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ ", marker, " + MinorityInd + HighRiskInd + risk_score"))
  
  fit <- survey::svycoxph(fm, design=design)
  
  dat <- tidy(fit, conf.int = T, exponentiate = T) %>% data.frame() %>%
    select(term, estimate, conf.low, conf.high, p.value) %>%
    filter(!term %in% c("MinorityInd", "HighRiskInd", "risk_score")) %>%
    rename(marker = term) %>%
    add_row(marker = "Lower", estimate = 1, conf.low = NA, conf.high = NA, p.value = NA) %>%
    mutate(marker_cut = c("Middle", "Upper", "Lower")) %>%
    arrange(marker_cut) %>%
    mutate(events = ifelse(marker_cut == "Lower", sum(data_lower$EventIndPrimaryD57),
                           ifelse(marker_cut == "Middle", sum(data_middle$EventIndPrimaryD57),
                                  sum(data_upper$EventIndPrimaryD57))),
           n = ifelse(marker_cut == "Lower", length(data_lower$EventIndPrimaryD57),
                      ifelse(marker_cut == "Middle", length(data_middle$EventIndPrimaryD57),
                             length(data_upper$EventIndPrimaryD57))),
           group = group)
  
}

################################################## 
# Day 57 analysis
dat.D57 <- dat %>% filter(!is.na(wt))

# Define trichotomized version of the markers
# trichotomize Day57bindSpike
vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D57$Day57bindSpike, weights = dat.D57$wt, probs = c(1/3, 2/3)), Inf)
if((vec_cutpoints[[2]] == min(dat.D57$Day57bindSpike, na.rm = T)) | 
   (vec_cutpoints[[3]] == max(dat.D57$Day57bindSpike, na.rm = T)) |
   (vec_cutpoints[[2]] == vec_cutpoints[[3]])){
  dat.D57$Day57bindSpikecat <- Hmisc::cut2(dat.D57$Day57bindSpike, breaks = 3)
} else {
  dat.D57$Day57bindSpikecat <- Hmisc::cut2(dat.D57$Day57bindSpike, vec_cutpoints)}
rm(vec_cutpoints)

# trichotomize Day57bindRBD
vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D57$Day57bindRBD, weights = dat.D57$wt, probs = c(1/3, 2/3)), Inf)
if((vec_cutpoints[[2]] == min(dat.D57$Day57bindRBD, na.rm = T)) | 
   (vec_cutpoints[[3]] == max(dat.D57$Day57bindRBD, na.rm = T)) |
   (vec_cutpoints[[2]] == vec_cutpoints[[3]])){
  dat.D57$Day57bindRBDcat <- Hmisc::cut2(dat.D57$Day57bindRBD, breaks = 3)
} else {
  dat.D57$Day57bindRBDcat <- Hmisc::cut2(dat.D57$Day57bindRBD, vec_cutpoints)}
rm(vec_cutpoints)

# trichotomize Day57pseudoneutid50
vec_cutpoints <- c(-Inf, Hmisc::wtd.quantile(dat.D57$Day57pseudoneutid50, weights = dat.D57$wt, probs = c(1/3, 2/3)), Inf)
if((vec_cutpoints[[2]] == min(dat.D57$Day57pseudoneutid50, na.rm = T)) | 
   (vec_cutpoints[[3]] == max(dat.D57$Day57pseudoneutid50, na.rm = T)) |
   (vec_cutpoints[[2]] == vec_cutpoints[[3]])){
  dat.D57$Day57pseudoneutid50cat <- Hmisc::cut2(dat.D57$Day57pseudoneutid50, breaks = 3)
} else {
  dat.D57$Day57pseudoneutid50cat <- Hmisc::cut2(dat.D57$Day57pseudoneutid50, vec_cutpoints)}
rm(vec_cutpoints)

# trichotomize Day57pseudoneutid80
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

##########################################


dat.D57.design <- twophase(id=list(~Ptid, ~Ptid),strata=list(NULL, ~Wstratum),
                           data=dat.D57, subset=~TwophasesampInd)

tab1 <- getHR_D57_continuous_marker(marker = "Day57bindSpike", design = dat.D57.design, data = dat.D57, group = "All baseline negative, vaccine") %>%
  bind_rows(getHR_D57_continuous_marker(marker = "Day57bindRBD", design = dat.D57.design, data = dat.D57, group = "All baseline negative, vaccine"),
            getHR_D57_continuous_marker(marker = "Day57pseudoneutid50", design = dat.D57.design, data = dat.D57, group = "All baseline negative, vaccine"),
            getHR_D57_continuous_marker(marker = "Day57pseudoneutid80", design = dat.D57.design, data = dat.D57, group = "All baseline negative, vaccine"))

tab2 <- getHR_D57_categorical_marker(marker = "Day57bindSpikecat", design = dat.D57.design, data = dat.D57, group = "All baseline negative, vaccine") %>%
  bind_rows(getHR_D57_categorical_marker(marker = "Day57bindRBDcat", design = dat.D57.design, data = dat.D57, group = "All baseline negative, vaccine"),
            getHR_D57_categorical_marker(marker = "Day57pseudoneutid50cat", design = dat.D57.design, data = dat.D57, group = "All baseline negative, vaccine"),
            getHR_D57_categorical_marker(marker = "Day57pseudoneutid80cat", design = dat.D57.design, data = dat.D57, group = "All baseline negative, vaccine")) 

################################################## 
# Age >= 65
# create data subsets
dat <- inputFile %>% filter(Trt == 1 & Bserostatus == 0 & Age >= 65)
dat.D57 <- dat %>% filter(!is.na(wt))




















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


