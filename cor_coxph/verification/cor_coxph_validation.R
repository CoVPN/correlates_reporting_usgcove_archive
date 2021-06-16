#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R")) 
#-----------------------------------------------

# load required libraries and functions
library(tidyverse)
library(here)
library(glmnet)
library(kyotil)
library(conflicted)
conflicted::conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
library(marginalizedRisk)
library(tools) 
library(survey)
library(Hmisc)
library(xtable) 
library(broom)

# Read in data file
inputFile <- read.csv(here::here("verification", "verification_input", "moderna_mock_data_processed_with_riskscore.csv")) %>%
  mutate(Day57bindSpike = ifelse(Day57bindSpike > log10(uloqs[["bindSpike"]]), log10(uloqs[["bindSpike"]]), Day57bindSpike), 
         Day57bindRBD = ifelse(Day57bindRBD > log10(uloqs[["bindRBD"]]), log10(uloqs[["bindRBD"]]), Day57bindRBD),
         Day57pseudoneutid50 = ifelse(Day57pseudoneutid50 > log10(uloqs[["pseudoneutid50"]]), log10(uloqs[["pseudoneutid50"]]), Day57pseudoneutid50),
         Day57pseudoneutid80 = ifelse(Day57pseudoneutid80 > log10(uloqs[["pseudoneutid80"]]), log10(uloqs[["pseudoneutid80"]]), Day57pseudoneutid80),
         
         Day29bindSpike = ifelse(Day29bindSpike > log10(uloqs[["bindSpike"]]), log10(uloqs[["bindSpike"]]), Day29bindSpike), 
         Day29bindRBD = ifelse(Day29bindRBD > log10(uloqs[["bindRBD"]]), log10(uloqs[["bindRBD"]]), Day29bindRBD),
         Day29pseudoneutid50 = ifelse(Day29pseudoneutid50 > log10(uloqs[["pseudoneutid50"]]), log10(uloqs[["pseudoneutid50"]]), Day29pseudoneutid50),
         Day29pseudoneutid80 = ifelse(Day29pseudoneutid80 > log10(uloqs[["pseudoneutid80"]]), log10(uloqs[["pseudoneutid80"]]), Day29pseudoneutid80)) 

# This function takes marker name as string, and the design.
# It returns the HR for the marker, CI and p.value
getHR_D57_continuous_marker <- function(marker, design, data, group){
  if(group %in% c("Age < 65, At risk", "Age < 65, Not at risk", "At risk", "Not at risk")){
    fm = as.formula(paste("Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ ", marker, " + MinorityInd + risk_score"))
  } else if (group %in% c("WhiteNonHispanic", "Comm. of color")){
    fm = as.formula(paste("Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ ", marker, " + HighRiskInd + risk_score"))
  }else{
    fm = as.formula(paste("Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ ", marker, " + MinorityInd + HighRiskInd + risk_score"))
  }
  
  fit <- survey::svycoxph(fm, design=design)
  
  tidy(fit, conf.int = T, exponentiate = T) %>% data.frame() %>%
    select(term, estimate, conf.low, conf.high, p.value) %>%
    filter(!term %in% c("MinorityInd", "HighRiskInd", "risk_score")) %>%
    rename(marker = term) %>%
    mutate(cases = sum(data$EventIndPrimaryD57),
           atRisk = length(data$EventIndPrimaryD57),
           group = group) 
}



getHR_D57_categorical_marker <- function(marker, design, data, group){
  
  # table(data %>% pull(all_of(marker)))
  # Add events and n for lower, middle and upper
  data_lower <- data %>% filter(get(marker) == (data %>% pull(marker) %>% levels())[1])
  data_middle <- data %>% filter(get(marker) == (data %>% pull(marker) %>% levels())[2])
  data_upper <- data %>% filter(get(marker) == (data %>% pull(marker) %>% levels())[3])
  
  fm = as.formula(paste("Surv(EventTimePrimaryD57, EventIndPrimaryD57) ~ ", marker, " + MinorityInd + HighRiskInd + risk_score"))
  
  fit <- survey::svycoxph(fm, design=design)
  
  dat <- tidy(fit, conf.int = T, exponentiate = T) %>% data.frame() %>%
    select(term, estimate, conf.low, conf.high, p.value) %>%
    filter(!term %in% c("MinorityInd", "HighRiskInd", "risk_score")) %>%
    rename(marker = term) %>%
    add_row(marker = "Lower", estimate = 1, conf.low = NA, conf.high = NA, p.value = NA) %>%
    mutate(marker_cut = c("Middle", "Upper", "Lower")) %>%
    arrange(marker_cut) %>%
    mutate(cases = ifelse(marker_cut == "Lower", sum(data_lower %>% filter(EventIndPrimaryD57 == 1) %>% .$wt.D57),
                           ifelse(marker_cut == "Middle", sum(data_middle %>% filter(EventIndPrimaryD57 == 1) %>% .$wt.D57),
                                  sum(data_upper %>% filter(EventIndPrimaryD57 == 1) %>% .$wt.D57))),
           atRisk = ifelse(marker_cut == "Lower", sum(data_lower %>% .$wt.D57),
                      ifelse(marker_cut == "Middle", sum(data_middle %>% .$wt.D57),
                             sum(data_upper %>% .$wt.D57))),
           `Attack rate` = cases/atRisk,
           group = group)
  
  # Generalized Wald test
  idx = 1:2
  stat=coef(fit)[idx] %*% solve(vcov(fit)[idx, idx]) %*% coef(fit)[idx] #  X' %*% V^-1 %*% X to get the chi-squared statistic
  dat <- dat %>% mutate(
    id = seq.int(nrow(.)), 
    overall.pval = ifelse(id == 1, format(round(pchisq(stat, length(idx), lower.tail = FALSE), 3), nsmall = 3), ''),
    overall.pval = ifelse(overall.pval < 0.001 & overall.pval!= '', "<0.001", overall.pval))
}

################################################## 
get_results_in_df_D57 <- function(dat, group){
  # Day 57 analysis
  dat.D57 <- dat %>% filter(!is.na(wt.D57))

  dat.D57.design <- twophase(id=list(~Ptid, ~Ptid),strata=list(NULL, ~Wstratum),
                             data=dat.D57, subset=~ph2.D57)
  
  tab1 <- getHR_D57_continuous_marker(marker = "Day57bindSpike", design = dat.D57.design, data = dat.D57, group = group) %>%
    bind_rows(getHR_D57_continuous_marker(marker = "Day57bindRBD", design = dat.D57.design, data = dat.D57, group = group),
              getHR_D57_continuous_marker(marker = "Day57pseudoneutid50", design = dat.D57.design, data = dat.D57, group = group),
              getHR_D57_continuous_marker(marker = "Day57pseudoneutid80", design = dat.D57.design, data = dat.D57, group = group)) %>%
    mutate(markerFormat = "continuous")
  
  if(group == "All baseline negative, vaccine"){
    
    # Define trichotomized version of the markers
    # trichotomize Day57bindSpike
    if(quantile(dat.D57$Day57bindSpike, probs = c(2/3), na.rm = T) > log10(uloqs[["bindSpike"]])*.9999){
      d = dat.D57 %>% filter(Day57bindSpike <= log10(uloqs[["bindSpike"]])*.9999)
      vec_cutpoints <- Hmisc::wtd.quantile(d$Day57bindSpike, weights = d$wt.D57, probs = c(1/2))
      vec_cutpoints[2] = log10(uloqs[["bindSpike"]])*.9999
    }else{
      vec_cutpoints <- Hmisc::wtd.quantile(dat.D57$Day57bindSpike, weights = dat.D57$wt.D57, probs = c(1/3, 2/3))
    }
    vec_cutpoints <- c(-Inf, vec_cutpoints, Inf)
    getCLASS = try(factor(cut(dat.D57$Day57bindSpike, breaks = vec_cutpoints)))
    if(class(getCLASS) == "factor"){
      dat.D57$Day57bindSpikecat = factor(cut(dat.D57$Day57bindSpike, breaks = vec_cutpoints))
      if(length(levels(dat.D57$Day57bindSpikecat))!=3){
        dat.D57$Day57bindSpikecat <- factor(cut(dat.D57$Day57bindSpike, breaks = 3))
      }
    }else if(class(getCLASS) == "try-error"){
      dat.D57$Day57bindSpikecat <- factor(cut(dat.D57$Day57bindSpike, breaks = 3))
    }
    rm(vec_cutpoints)
    
    # trichotomize Day57bindRBD
    if(quantile(dat.D57$Day57bindRBD, probs = c(2/3), na.rm = T) > log10(uloqs[["bindRBD"]])*.9999){
      d = dat.D57 %>% filter(Day57bindRBD <= log10(uloqs[["bindRBD"]])*.9999)
      vec_cutpoints <- Hmisc::wtd.quantile(d$Day57bindRBD, weights = d$wt.D57, probs = c(1/2))
      vec_cutpoints[2] = log10(uloqs[["bindRBD"]])*.9999
    }else{
      vec_cutpoints <- Hmisc::wtd.quantile(dat.D57$Day57bindRBD, weights = dat.D57$wt.D57, probs = c(1/3, 2/3))
    }
    vec_cutpoints <- c(-Inf, vec_cutpoints, Inf)
    getCLASS = try(factor(cut(dat.D57$Day57bindRBD, breaks = vec_cutpoints)))
    if(class(getCLASS) == "factor"){
      dat.D57$Day57bindRBDcat = factor(cut(dat.D57$Day57bindRBD, breaks = vec_cutpoints))
      if(length(levels(dat.D57$Day57bindRBDcat))!=3){
        dat.D57$Day57bindRBDcat <- factor(cut(dat.D57$Day57bindRBD, breaks = 3))
      }
    }else if(class(getCLASS) == "try-error"){
      dat.D57$Day57bindRBDcat <- factor(cut(dat.D57$Day57bindRBD, breaks = 3))
    }
    rm(vec_cutpoints)
    
    # trichotomize Day57pseudoneutid50
    if(quantile(dat.D57$Day57pseudoneutid50, probs = c(2/3), na.rm = T) > log10(uloqs[["pseudoneutid50"]])*.9999){
      d = dat.D57 %>% filter(Day57pseudoneutid50 <= log10(uloqs[["pseudoneutid50"]])*.9999)
      vec_cutpoints <- Hmisc::wtd.quantile(d$Day57pseudoneutid50, weights = d$wt.D57, probs = c(1/2))
      vec_cutpoints[2] = log10(uloqs[["pseudoneutid50"]])*.9999
    }else{
      vec_cutpoints <- Hmisc::wtd.quantile(dat.D57$Day57pseudoneutid50, weights = dat.D57$wt.D57, probs = c(1/3, 2/3))
    }
    vec_cutpoints <- c(-Inf, vec_cutpoints, Inf)
    getCLASS = try(factor(cut(dat.D57$Day57pseudoneutid50, breaks = vec_cutpoints)))
    if(class(getCLASS) == "factor"){
      dat.D57$Day57pseudoneutid50cat = factor(cut(dat.D57$Day57pseudoneutid50, breaks = vec_cutpoints))
      if(length(levels(dat.D57$Day57pseudoneutid50cat))!=3){
        dat.D57$Day57pseudoneutid50cat <- factor(cut(dat.D57$Day57pseudoneutid50, breaks = 3))
      }
    }else if(class(getCLASS) == "try-error"){
      dat.D57$Day57pseudoneutid50cat <- factor(cut(dat.D57$Day57pseudoneutid50, breaks = 3))
    }
    rm(vec_cutpoints)
    
    # trichotomize Day57pseudoneutid80
    if(quantile(dat.D57$Day57pseudoneutid80, probs = c(2/3), na.rm = T) > log10(uloqs[["pseudoneutid80"]])*.9999){
      d = dat.D57 %>% filter(Day57pseudoneutid80 <= log10(uloqs[["pseudoneutid80"]])*.9999)
      vec_cutpoints <- Hmisc::wtd.quantile(d$Day57pseudoneutid80, weights = d$wt.D57, probs = c(1/2))
      vec_cutpoints[2] = log10(uloqs[["pseudoneutid80"]])*.9999
    }else{
      vec_cutpoints <- Hmisc::wtd.quantile(dat.D57$Day57pseudoneutid80, weights = dat.D57$wt.D57, probs = c(1/3, 2/3))
    }
    vec_cutpoints <- c(-Inf, vec_cutpoints, Inf)
    getCLASS = try(factor(cut(dat.D57$Day57pseudoneutid80, breaks = vec_cutpoints)))
    if(class(getCLASS) == "factor"){
      dat.D57$Day57pseudoneutid80cat = factor(cut(dat.D57$Day57pseudoneutid80, breaks = vec_cutpoints))
      if(length(levels(dat.D57$Day57pseudoneutid80cat))!=3){
        dat.D57$Day57pseudoneutid80cat <- factor(cut(dat.D57$Day57pseudoneutid80, breaks = 3))
      }
    }else if(class(getCLASS) == "try-error"){
      dat.D57$Day57pseudoneutid80cat <- factor(cut(dat.D57$Day57pseudoneutid80, breaks = 3))
    }
    rm(vec_cutpoints)
    
    dat.D57.design <- twophase(id=list(~Ptid, ~Ptid),strata=list(NULL, ~Wstratum),
                               data=dat.D57, subset=~ph2.D57)
    
    tab2 <- getHR_D57_categorical_marker(marker = "Day57bindSpikecat", design = dat.D57.design, data = dat.D57, group = group) %>%
      bind_rows(getHR_D57_categorical_marker(marker = "Day57bindRBDcat", design = dat.D57.design, data = dat.D57, group = group),
                getHR_D57_categorical_marker(marker = "Day57pseudoneutid50cat", design = dat.D57.design, data = dat.D57, group = group),
                getHR_D57_categorical_marker(marker = "Day57pseudoneutid80cat", design = dat.D57.design, data = dat.D57, group = group)) %>%
      mutate(markerFormat = "categorical")
    
    return(bind_rows(tab1, tab2))
  }
  
  if(group != "All baseline negative, vaccine"){
    return(tab1)
  }
}


################################################## 
# All D57 analyses
tab <- get_results_in_df_D57(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0), group = "All baseline negative, vaccine") %>%
  bind_rows(get_results_in_df_D57(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & Age >= 65), group = "Age >= 65"),
            get_results_in_df_D57(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & Age < 65 & HighRiskInd == 1), group = "Age < 65, At risk"),
            get_results_in_df_D57(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & Age < 65 & HighRiskInd == 0), group = "Age < 65, Not at risk"),
            get_results_in_df_D57(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & Age < 65), group = "Age < 65"),
            get_results_in_df_D57(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & HighRiskInd == 1), group = "At risk"),
            get_results_in_df_D57(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & HighRiskInd == 0), group = "Not at risk"),
            get_results_in_df_D57(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & MinorityInd==1), group = "Comm. of color"),
            get_results_in_df_D57(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & MinorityInd == 0), group = "WhiteNonHispanic"),
            get_results_in_df_D57(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & Sex == 1), group = "Men"),
            get_results_in_df_D57(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & Sex == 0), group = "Women"))

tab %>% 
  write.csv(here("verification", "verification_output", "D57.mock.csv"))

tab <- tab %>% select(marker, cases, atRisk, estimate, conf.low, conf.high, p.value)

###############################################################################
###############################################################################
###############################################################################

# This function takes marker name as string, and the design.
# It returns the HR for the marker, CI and p.value
getHR_D29_continuous_marker <- function(marker, design, data, group){
  if(group %in% c("Age < 65, At risk", "Age < 65, Not at risk", "At risk", "Not at risk")){
    fm = as.formula(paste("Surv(EventTimePrimaryD29, EventIndPrimaryD29) ~ ", marker, " + MinorityInd + risk_score"))
  } else if (group %in% c("WhiteNonHispanic", "Comm. of color")){
    fm = as.formula(paste("Surv(EventTimePrimaryD29, EventIndPrimaryD29) ~ ", marker, " + HighRiskInd + risk_score"))
  }else{
    fm = as.formula(paste("Surv(EventTimePrimaryD29, EventIndPrimaryD29) ~ ", marker, " + MinorityInd + HighRiskInd + risk_score"))
  }
  
  fit <- survey::svycoxph(fm, design=design)
  
  tidy(fit, conf.int = T, exponentiate = T) %>% data.frame() %>%
    select(term, estimate, conf.low, conf.high, p.value) %>%
    filter(!term %in% c("MinorityInd", "HighRiskInd", "risk_score")) %>%
    rename(marker = term) %>%
    mutate(cases = sum(data$EventIndPrimaryD29),
           atRisk = length(data$EventIndPrimaryD29),
           group = group) 
}



getHR_D29_categorical_marker <- function(marker, design, data, group){
  
  # table(data %>% pull(all_of(marker)))
  # Add events and n for lower, middle and upper
  data_lower <- data %>% filter(get(marker) == (data %>% pull(marker) %>% levels())[1])
  data_middle <- data %>% filter(get(marker) == (data %>% pull(marker) %>% levels())[2])
  data_upper <- data %>% filter(get(marker) == (data %>% pull(marker) %>% levels())[3])
  
  fm = as.formula(paste("Surv(EventTimePrimaryD29, EventIndPrimaryD29) ~ ", marker, " + MinorityInd + HighRiskInd + risk_score"))
  
  fit <- survey::svycoxph(fm, design=design)
  
  dat <- tidy(fit, conf.int = T, exponentiate = T) %>% data.frame() %>%
    select(term, estimate, conf.low, conf.high, p.value) %>%
    filter(!term %in% c("MinorityInd", "HighRiskInd", "risk_score")) %>%
    rename(marker = term) %>%
    add_row(marker = "Lower", estimate = 1, conf.low = NA, conf.high = NA, p.value = NA) %>%
    mutate(marker_cut = c("Middle", "Upper", "Lower")) %>%
    arrange(marker_cut) %>%
    mutate(cases = ifelse(marker_cut == "Lower", sum(data_lower %>% filter(EventIndPrimaryD29 == 1) %>% .$wt.D29),
                          ifelse(marker_cut == "Middle", sum(data_middle %>% filter(EventIndPrimaryD29 == 1) %>% .$wt.D29),
                                 sum(data_upper %>% filter(EventIndPrimaryD29 == 1) %>% .$wt.D29))),
           atRisk = ifelse(marker_cut == "Lower", sum(data_lower %>% .$wt.D29),
                           ifelse(marker_cut == "Middle", sum(data_middle %>% .$wt.D29),
                                  sum(data_upper %>% .$wt.D29))),
           `Attack rate` = round(cases, 0)/atRisk,
           group = group)
  
  # Generalized Wald test
  idx = 1:2
  stat=coef(fit)[idx] %*% solve(vcov(fit)[idx, idx]) %*% coef(fit)[idx] #  X' %*% V^-1 %*% X to get the chi-squared statistic
  dat <- dat %>% mutate(
    id = seq.int(nrow(.)), 
    overall.pval = ifelse(id == 1, format(round(pchisq(stat, length(idx), lower.tail = FALSE), 3), nsmall = 3), ''),
    overall.pval = ifelse(overall.pval < 0.001 & overall.pval!= '', "<0.001", overall.pval))
}

################################################## 
get_results_in_df_D29 <- function(dat, group){
  # Day 29 analysis 
  dat.D29 <- dat %>% filter(!is.na(wt.D29))
  
  dat.D29.design <- twophase(id=list(~Ptid, ~Ptid),strata=list(NULL, ~Wstratum),
                             data=dat.D29, subset=~ph2.D29)
  
  tab1 <- getHR_D29_continuous_marker(marker = "Day29bindSpike", design = dat.D29.design, data = dat.D29, group = group) %>%
    bind_rows(getHR_D29_continuous_marker(marker = "Day29bindRBD", design = dat.D29.design, data = dat.D29, group = group),
              getHR_D29_continuous_marker(marker = "Day29pseudoneutid50", design = dat.D29.design, data = dat.D29, group = group),
              getHR_D29_continuous_marker(marker = "Day29pseudoneutid80", design = dat.D29.design, data = dat.D29, group = group)) %>%
    mutate(markerFormat = "continuous")
  
  if(group == "All baseline negative, vaccine"){
    
    # Define trichotomized version of the markers
    # trichotomize Day29bindSpike
    if(quantile(dat.D29$Day29bindSpike, probs = c(2/3), na.rm = T) > log10(uloqs[["bindSpike"]])*.9999){
      d = dat.D29 %>% filter(Day29bindSpike <= log10(uloqs[["bindSpike"]])*.9999)
      vec_cutpoints <- Hmisc::wtd.quantile(d$Day29bindSpike, weights = d$wt.D29, probs = c(1/2))
      vec_cutpoints[2] = log10(uloqs[["bindSpike"]])*.9999
    }else{
      vec_cutpoints <- Hmisc::wtd.quantile(dat.D29$Day29bindSpike, weights = dat.D29$wt.D29, probs = c(1/3, 2/3))
    }
    vec_cutpoints <- c(-Inf, vec_cutpoints, Inf)
    getCLASS = try(factor(cut(dat.D29$Day29bindSpike, breaks = vec_cutpoints)))
    if(class(getCLASS) == "factor"){
      dat.D29$Day29bindSpikecat = factor(cut(dat.D29$Day29bindSpike, breaks = vec_cutpoints))
      if(length(levels(dat.D29$Day29bindSpikecat))!=3){
        dat.D29$Day29bindSpikecat <- factor(cut(dat.D29$Day29bindSpike, breaks = 3))
      }
    }else if(class(getCLASS) == "try-error"){
      dat.D29$Day29bindSpikecat <- factor(cut(dat.D29$Day29bindSpike, breaks = 3))
    }
    rm(vec_cutpoints)
    
    # trichotomize Day29bindRBD
    if(quantile(dat.D29$Day29bindRBD, probs = c(2/3), na.rm = T) > log10(uloqs[["bindRBD"]])*.9999){
      d = dat.D29 %>% filter(Day29bindRBD <= log10(uloqs[["bindRBD"]])*.9999)
      vec_cutpoints <- Hmisc::wtd.quantile(d$Day29bindRBD, weights = d$wt.D29, probs = c(1/2))
      vec_cutpoints[2] = log10(uloqs[["bindRBD"]])*.9999
    }else{
      vec_cutpoints <- Hmisc::wtd.quantile(dat.D29$Day29bindRBD, weights = dat.D29$wt.D29, probs = c(1/3, 2/3))
    }
    vec_cutpoints <- c(-Inf, vec_cutpoints, Inf)
    getCLASS = try(factor(cut(dat.D29$Day29bindRBD, breaks = vec_cutpoints)))
    if(class(getCLASS) == "factor"){
      dat.D29$Day29bindRBDcat = factor(cut(dat.D29$Day29bindRBD, breaks = vec_cutpoints))
      if(length(levels(dat.D29$Day29bindRBDcat))!=3){
        dat.D29$Day29bindRBDcat <- factor(cut(dat.D29$Day29bindRBD, breaks = 3))
      }
    }else if(class(getCLASS) == "try-error"){
      dat.D29$Day29bindRBDcat <- factor(cut(dat.D29$Day29bindRBD, breaks = 3))
    }
    rm(vec_cutpoints)
    
    # trichotomize Day29pseudoneutid50
    if(quantile(dat.D29$Day29pseudoneutid50, probs = c(2/3), na.rm = T) > log10(uloqs[["pseudoneutid50"]])*.9999){
      d = dat.D29 %>% filter(Day29pseudoneutid50 <= log10(uloqs[["pseudoneutid50"]])*.9999)
      vec_cutpoints <- Hmisc::wtd.quantile(d$Day29pseudoneutid50, weights = d$wt.D29, probs = c(1/2))
      vec_cutpoints[2] = log10(uloqs[["pseudoneutid50"]])*.9999
    }else{
      vec_cutpoints <- Hmisc::wtd.quantile(dat.D29$Day29pseudoneutid50, weights = dat.D29$wt.D29, probs = c(1/3, 2/3))
    }
    vec_cutpoints <- c(-Inf, vec_cutpoints, Inf)
    getCLASS = try(factor(cut(dat.D29$Day29pseudoneutid50, breaks = vec_cutpoints)))
    if(class(getCLASS) == "factor"){
      dat.D29$Day29pseudoneutid50cat = factor(cut(dat.D29$Day29pseudoneutid50, breaks = vec_cutpoints))
      if(length(levels(dat.D29$Day29pseudoneutid50cat))!=3){
        dat.D29$Day29pseudoneutid50cat <- factor(cut(dat.D29$Day29pseudoneutid50, breaks = 3))
      }
    }else if(class(getCLASS) == "try-error"){
      dat.D29$Day29pseudoneutid50cat <- factor(cut(dat.D29$Day29pseudoneutid50, breaks = 3))
    }
    rm(vec_cutpoints)
    
    # trichotomize Day29pseudoneutid80
    if(quantile(dat.D29$Day29pseudoneutid80, probs = c(2/3), na.rm = T) > log10(uloqs[["pseudoneutid80"]])*.9999){
      d = dat.D29 %>% filter(Day29pseudoneutid80 <= log10(uloqs[["pseudoneutid80"]])*.9999)
      vec_cutpoints <- Hmisc::wtd.quantile(d$Day29pseudoneutid80, weights = d$wt.D29, probs = c(1/2))
      vec_cutpoints[2] = log10(uloqs[["pseudoneutid80"]])*.9999
    }else{
      vec_cutpoints <- Hmisc::wtd.quantile(dat.D29$Day29pseudoneutid80, weights = dat.D29$wt.D29, probs = c(1/3, 2/3))
    }
    vec_cutpoints <- c(-Inf, vec_cutpoints, Inf)
    getCLASS = try(factor(cut(dat.D29$Day29pseudoneutid80, breaks = vec_cutpoints)))
    if(class(getCLASS) == "factor"){
      dat.D29$Day29pseudoneutid80cat = factor(cut(dat.D29$Day29pseudoneutid80, breaks = vec_cutpoints))
      if(length(levels(dat.D29$Day29pseudoneutid80cat))!=3){
        dat.D29$Day29pseudoneutid80cat <- factor(cut(dat.D29$Day29pseudoneutid80, breaks = 3))
      }
    }else if(class(getCLASS) == "try-error"){
      dat.D29$Day29pseudoneutid80cat <- factor(cut(dat.D29$Day29pseudoneutid80, breaks = 3))
    }
    rm(vec_cutpoints)
    
    dat.D29.design <- twophase(id=list(~Ptid, ~Ptid),strata=list(NULL, ~Wstratum),
                               data=dat.D29, subset=~ph2.D29)
    
    tab2 <- getHR_D29_categorical_marker(marker = "Day29bindSpikecat", design = dat.D29.design, data = dat.D29, group = group) %>%
      bind_rows(getHR_D29_categorical_marker(marker = "Day29bindRBDcat", design = dat.D29.design, data = dat.D29, group = group),
                getHR_D29_categorical_marker(marker = "Day29pseudoneutid50cat", design = dat.D29.design, data = dat.D29, group = group),
                getHR_D29_categorical_marker(marker = "Day29pseudoneutid80cat", design = dat.D29.design, data = dat.D29, group = group)) %>%
      mutate(markerFormat = "categorical")
    
    return(bind_rows(tab1, tab2))
  }
  
  if(group != "All baseline negative, vaccine"){
    return(tab1)
  }
}


################################################## 
# All D29 analyses
tab <- get_results_in_df_D29(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0), group = "All baseline negative, vaccine") %>%
  bind_rows(get_results_in_df_D29(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & Age >= 65), group = "Age >= 65"),
            get_results_in_df_D29(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & Age < 65 & HighRiskInd == 1), group = "Age < 65, At risk"),
            get_results_in_df_D29(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & Age < 65 & HighRiskInd == 0), group = "Age < 65, Not at risk"),
            get_results_in_df_D29(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & Age < 65), group = "Age < 65"),
            get_results_in_df_D29(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & HighRiskInd == 1), group = "At risk"),
            get_results_in_df_D29(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & HighRiskInd == 0), group = "Not at risk"),
            get_results_in_df_D29(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & MinorityInd==1), group = "Comm. of color"),
            get_results_in_df_D29(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & MinorityInd == 0), group = "WhiteNonHispanic"),
            get_results_in_df_D29(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & Sex == 1), group = "Men"),
            get_results_in_df_D29(dat = inputFile %>% filter(Trt == 1 & Bserostatus == 0 & Sex == 0), group = "Women"))

tab %>% 
  write.csv(here("verification", "verification_output", "D29.mock.csv"))
